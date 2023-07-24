const std = @import("std");
const ray = @cImport({
    @cInclude("raylib.h");
});

var prng = std.rand.DefaultPrng.init(42);
var rng = prng.random();

const Vec2 = @Vector(2, f32);
const Vec2i = @Vector(2, c_int);

const SCREEN_WIDTH = 800.0;
const SCREEN_HEIGHT = 600.0;

const UP = SCREEN_HEIGHT / SCREEN_WIDTH;
const DOWN = -SCREEN_HEIGHT / SCREEN_WIDTH;

const Balls = struct {
    pos: [*]Vec2,
    vel: [*]Vec2,
    force: [*]Vec2,
    mass: [*]f32,
    radius: [*]f32,
    color: [*]ray.Color,
    count: usize,
};

fn norm_to_pixel_v2(pos: Vec2) ray.Vector2 {
    return ray.Vector2{
        .x = (pos[0] + 1) * SCREEN_WIDTH / 2,
        .y = (-pos[1] * SCREEN_WIDTH + SCREEN_HEIGHT) / 2,
    };
}

fn norm_to_pixel_f(pos: f32) f32 {
    return pos * SCREEN_WIDTH / 2;
}

fn dist_v2(v0: Vec2, v1: Vec2) f32 {
    var dx = v0[0] - v1[0];
    var dy = v0[1] - v1[1];
    return abs_v2(Vec2{ dx, dy });
}

fn dist_square_v2(v0: Vec2, v1: Vec2) f32 {
    var dx = v0[0] - v1[0];
    var dy = v0[1] - v1[1];
    return abs_square_v2(Vec2{ dx, dy });
}

fn abs_f32(v: f32) f32 {
    return if (v < 0) -v else v;
}

fn abs_v2(v: Vec2) f32 {
    var x = v[0];
    var y = v[1];
    var m = std.math.max(abs_f32(x), abs_f32(y));
    if (m == 0) {
        return 0;
    } else {
        x /= m;
        y /= m;
        return m * std.math.sqrt(x * x + y * y);
    }
}

fn abs_square_v2(v: Vec2) f32 {
    return v[0] * v[0] + v[1] * v[1];
}

fn rand_range(min: f32, max: f32) f32 {
    return @intToFloat(f32, rng.int(u16)) / @as(f32, std.math.maxInt(u16)) * (max - min) + min;
}

fn normalize_v2(v: *Vec2) void {
    var abs = abs_v2(v.*);
    v.* /= @splat(2, abs);
}

fn cross(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[1] - v0[1] * v1[0];
}

fn dot(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[0] + v0[1] * v1[1];
}

fn intersect_ball_against_walls(balls: *Balls, index: usize) void {
    var acc = Vec2{ 0, 0 };
    var pos = balls.pos[index];
    var radius = balls.radius[index];

    if (pos[0] - radius < -1) {
        acc[0] -= pos[0] - radius + 1;
        balls.force[index][0] = 0;
    }

    if (pos[0] + radius > 1) {
        acc[0] -= pos[0] + radius - 1;
        balls.force[index][0] = 0;
    }

    if (pos[1] - radius < DOWN) {
        acc[1] -= pos[1] - radius - DOWN;
        balls.force[index][1] = 0;
    }

    if (pos[1] + radius > UP) {
        acc[1] -= pos[1] + radius - UP;
        balls.force[index][1] = 0;
    }

    balls.pos[index] += acc;
}

pub fn main() void {
    ray.InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello, Traingles!");
    ray.SetTargetFPS(60);

    var _balls_pos: [3]Vec2 = undefined;
    var _balls_vel: [3]Vec2 = undefined;
    var _balls_force: [3]Vec2 = undefined;
    var _balls_mass = [3]f32{ 2, 7, 10 };
    var _balls_color = [3]ray.Color{ ray.MAGENTA, ray.BLUE, ray.DARKBLUE };
    var _balls_radius = [3]f32{ 0.1, 0.2, 0.07 };

    var balls = Balls{
        .pos = &_balls_pos,
        .vel = &_balls_vel,
        .force = &_balls_force,
        .mass = &_balls_mass,
        .radius = &_balls_radius,
        .color = &_balls_color,
        .count = 3,
    };

    {
        var i: usize = 0;
        while (i < balls.count) : (i += 1) {
            balls.pos[i][0] = rand_range(-0.9, 0.9);
            balls.pos[i][1] = rand_range(DOWN + 0.1, UP - 0.1);

            balls.vel[i][0] = 0;
            balls.vel[i][1] = 0;
        }
    }

    const dt: f32 = 1.0 / 60.0;

    while (!ray.WindowShouldClose()) {
        ray.BeginDrawing();
        {
            const GRAVITY_CONSTANT = 0.005;

            var i: usize = 0;
            while (i < balls.count) : (i += 1) {
                balls.force[i][0] = 0;
                balls.force[i][1] = 0;
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                var j: usize = i + 1;
                while (j < balls.count) : (j += 1) {
                    var disp = balls.pos[j] - balls.pos[i];
                    var disp_len = abs_v2(disp);
                    var force = GRAVITY_CONSTANT * balls.mass[i] * balls.mass[j] / (disp_len * disp_len);

                    disp /= @splat(2, disp_len);

                    balls.force[i] += disp * @splat(2, force);
                    balls.force[j] -= disp * @splat(2, force);
                }
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                intersect_ball_against_walls(&balls, i);

                var j: usize = i + 1;
                while (j < balls.count) : (j += 1) {
                    var disp = balls.pos[j] - balls.pos[i];
                    var disp_len = abs_v2(disp);

                    disp /= @splat(2, disp_len);

                    var overlap = balls.radius[i] + balls.radius[j] - disp_len;
                    if (overlap > 0) {
                        overlap /= 2;
                        balls.pos[i] -= @splat(2, overlap) * disp;
                        balls.pos[j] += @splat(2, overlap) * disp;

                        var perp = if (cross(balls.force[i], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls.force[i] = perp * @splat(2, dot(balls.force[i], perp));

                        perp = if (cross(balls.force[j], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls.force[j] = perp * @splat(2, dot(balls.force[j], perp));
                    }
                }
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                balls.vel[i] += balls.force[i] / @splat(2, balls.mass[i]) * @splat(2, dt);
                balls.pos[i] += balls.vel[i] * @splat(2, dt);

                var force = balls.force[i];
                normalize_v2(&force);
                force *= @splat(2, @as(f32, 0.2));
                force += balls.pos[i];

                var center = norm_to_pixel_v2(balls.pos[i]);
                ray.DrawCircleV(center, norm_to_pixel_f(balls.radius[i]), balls.color[i]);
                ray.DrawLineEx(center, norm_to_pixel_v2(force), 3, ray.RED);
            }

            ray.ClearBackground(ray.DARKPURPLE);
        }
        ray.EndDrawing();
    }

    ray.CloseWindow();
}
