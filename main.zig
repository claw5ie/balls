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
    disp: [*]Vec2,
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

fn normalize_v2(v: Vec2) Vec2 {
    var result = v;
    var len = abs_v2(v);
    result /= @splat(2, len);

    return result;
}

fn cross(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[1] - v0[1] * v1[0];
}

fn dot(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[0] + v0[1] * v1[1];
}

fn reflect(vec: Vec2, dir: Vec2) Vec2 {
    var d = normalize_v2(dir);
    var w = Vec2{ cross(vec, d), dot(d, vec) };

    return Vec2{ cross(d, w), dot(d, w) };
}

fn intersect_ball_against_walls(balls: *Balls, index: usize) void {
    var pos = balls.pos[index];
    var radius = balls.radius[index];

    var overlap = pos[0] - radius + 1;
    if (overlap < 0) {
        balls.disp[index][0] -= overlap;
        if (balls.force[index][0] < 0) {
            balls.force[index][0] = 0;
        }

        if (balls.vel[index][0] < 0) {
            balls.vel[index][0] = -balls.vel[index][0];
        }
    }

    overlap = pos[0] + radius - 1;
    if (overlap > 0) {
        balls.disp[index][0] -= overlap;
        if (balls.force[index][0] > 0) {
            balls.force[index][0] = 0;
        }

        if (balls.vel[index][0] > 0) {
            balls.vel[index][0] = -balls.vel[index][0];
        }
    }

    overlap = pos[1] - radius - DOWN;
    if (overlap < 0) {
        balls.disp[index][1] -= overlap;
        if (balls.force[index][1] < 0) {
            balls.force[index][1] = 0;
        }

        if (balls.vel[index][1] < 0) {
            balls.vel[index][1] = -balls.vel[index][1];
        }
    }

    overlap = pos[1] + radius - UP;
    if (overlap > 0) {
        balls.disp[index][1] -= overlap;
        if (balls.force[index][1] > 0) {
            balls.force[index][1] = 0;
        }

        if (balls.vel[index][1] > 0) {
            balls.vel[index][1] = -balls.vel[index][1];
        }
    }
}

fn add_gravity_force(balls: *Balls, i: usize, j: usize) void {
    const GRAVITY_CONSTANT = 0.005;

    var disp = balls.pos[j] - balls.pos[i];
    var disp_len = abs_v2(disp);
    var force = GRAVITY_CONSTANT * balls.mass[i] * balls.mass[j] / (disp_len * disp_len);

    disp /= @splat(2, disp_len);

    balls.force[i] += disp * @splat(2, force);
    balls.force[j] -= disp * @splat(2, force);
}

fn add_spring_force(balls: *Balls, i: usize, j: usize) void {
    const STIFFNESS = 1.0;
    const REST_LENGTH = 0.1;
    const DAMPING_FACTOR = 0.1;

    var disp = balls.pos[j] - balls.pos[i];
    var disp_len = abs_v2(disp);

    disp /= @splat(2, disp_len);

    var damping = dot(disp, balls.vel[j] - balls.vel[i]) * DAMPING_FACTOR;
    var force = STIFFNESS * (disp_len - REST_LENGTH) + damping;

    balls.force[i] += @splat(2, force) * disp;
    balls.force[j] -= @splat(2, force) * disp;
}

pub fn main() void {
    ray.InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello, Traingles!");
    ray.SetTargetFPS(60);

    const MAX_BALL_COUNT = 3;
    var _balls_pos: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_vel: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_force: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_disp: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_mass = [MAX_BALL_COUNT]f32{ 2, 7, 5 };
    var _balls_color = [MAX_BALL_COUNT]ray.Color{
        ray.MAGENTA,
        ray.BLUE,
        ray.DARKBLUE,
    };
    var _balls_radius = [MAX_BALL_COUNT]f32{ 0.05, 0.05, 0.15 };

    var balls = Balls{
        .pos = &_balls_pos,
        .vel = &_balls_vel,
        .force = &_balls_force,
        .disp = &_balls_disp,
        .mass = &_balls_mass,
        .radius = &_balls_radius,
        .color = &_balls_color,
        .count = 3,
    };

    {
        var i: usize = 0;
        while (i < balls.count) : (i += 1) {
            balls.pos[i] = .{
                rand_range(-0.9, 0.9),
                rand_range(DOWN + 0.1, UP - 0.1),
            };
            balls.vel[i] = .{ 0, 0 };
        }
    }

    const dt: f32 = 1.0 / 60.0;

    while (!ray.WindowShouldClose()) {
        ray.BeginDrawing();
        {
            var i: usize = 0;
            while (i < balls.count) : (i += 1) {
                balls.force[i] = .{ 0, 0 };
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                var j: usize = i + 1;
                while (j < balls.count) : (j += 1) {
                    add_spring_force(&balls, i, j);
                }
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                balls.disp[i] = .{ 0, 0 };
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

                        balls.disp[i] -= @splat(2, overlap) * disp;
                        balls.disp[j] += @splat(2, overlap) * disp;

                        var perp = if (cross(balls.force[i], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls.force[i] = perp * @splat(2, dot(balls.force[i], perp));
                        if (dot(balls.vel[i], disp) > 0) {
                            balls.vel[i] = reflect(balls.vel[i], perp);
                        }

                        perp = if (cross(balls.force[j], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls.force[j] = perp * @splat(2, dot(balls.force[j], perp));
                        if (dot(balls.vel[j], disp) < 0) {
                            balls.vel[j] = reflect(balls.vel[j], perp);
                        }
                    }
                }
            }

            i = 0;
            while (i < balls.count) : (i += 1) {
                balls.pos[i] += balls.disp[i];

                balls.vel[i] += balls.force[i] / @splat(2, balls.mass[i]) * @splat(2, dt);
                balls.pos[i] += balls.vel[i] * @splat(2, dt);

                var center = norm_to_pixel_v2(balls.pos[i]);

                ray.DrawCircleV(center, norm_to_pixel_f(balls.radius[i]), balls.color[i]);
                ray.DrawLineEx(center, norm_to_pixel_v2(balls.pos[i] + balls.force[i] / @splat(2, @as(f32, 4))), 3, ray.RED);
                ray.DrawLineEx(center, norm_to_pixel_v2(balls.pos[i] + balls.vel[i] / @splat(2, @as(f32, 4))), 3, ray.GREEN);
            }

            ray.ClearBackground(ray.DARKPURPLE);
        }
        ray.EndDrawing();
    }

    ray.CloseWindow();
}
