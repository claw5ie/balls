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

const up = SCREEN_HEIGHT / SCREEN_WIDTH / 2.0;
const down = -SCREEN_HEIGHT / SCREEN_WIDTH / 2.0;

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

pub fn main() void {
    ray.InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello, Traingles!");
    ray.SetTargetFPS(60);

    const BALL_COUNT = 3;
    var balls_pos: [BALL_COUNT]Vec2 = undefined;
    var balls_vel: [BALL_COUNT]Vec2 = undefined;
    var balls_force: [BALL_COUNT]Vec2 = undefined;
    const balls_mass = [BALL_COUNT]f32{ 1, 2, 30 };
    const balls_color = [BALL_COUNT]ray.Color{ ray.MAGENTA, ray.BLUE, ray.DARKBLUE };
    var balls_radius = [BALL_COUNT]f32{ 0.1, 0.2, 0.02 };

    for (balls_pos) |*pos| {
        pos.*[0] = rand_range(-0.9, 0.9);
        pos.*[1] = rand_range(down + 0.1, up - 0.1);
    }

    std.mem.set(Vec2, &balls_vel, Vec2{ 0, 0 });

    const dt: f32 = 1.0 / 60.0;

    while (!ray.WindowShouldClose()) {
        ray.BeginDrawing();
        {
            const GRAVITY_CONSTANT = 0.001;

            std.mem.set(Vec2, &balls_force, Vec2{ 0, 0 });
            var i: usize = 0;
            while (i < BALL_COUNT) : (i += 1) {
                var j: usize = i + 1;
                while (j < BALL_COUNT) : (j += 1) {
                    var disp = balls_pos[j] - balls_pos[i];
                    var disp_len = abs_v2(disp);
                    var force = GRAVITY_CONSTANT * balls_mass[i] * balls_mass[j] / (disp_len * disp_len);

                    disp /= @splat(2, disp_len);

                    balls_force[i] += disp * @splat(2, force);
                    balls_force[j] -= disp * @splat(2, force);
                }
            }

            i = 0;
            while (i < BALL_COUNT) : (i += 1) {
                var j: usize = i + 1;
                while (j < BALL_COUNT) : (j += 1) {
                    var disp = balls_pos[j] - balls_pos[i];
                    var disp_len = abs_v2(disp);

                    disp /= @splat(2, disp_len);

                    var half_overlap = balls_radius[i] + balls_radius[j];
                    if (half_overlap > disp_len) {
                        half_overlap = (half_overlap - disp_len) / 2;
                        balls_pos[i] -= @splat(2, half_overlap) * disp;
                        balls_pos[j] += @splat(2, half_overlap) * disp;

                        var perp = if (cross(balls_force[i], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls_force[i] = perp * @splat(2, dot(balls_force[i], perp));

                        perp = if (cross(balls_force[j], disp) > 0)
                            Vec2{ disp[1], -disp[0] }
                        else
                            Vec2{ -disp[1], disp[0] };

                        balls_force[j] = perp * @splat(2, dot(balls_force[j], perp));
                    }
                }
            }

            i = 0;
            while (i < BALL_COUNT) : (i += 1) {
                balls_vel[i] += balls_force[i] / @splat(2, balls_mass[i]) * @splat(2, dt);
                balls_pos[i] += balls_vel[i] * @splat(2, dt);

                var force = balls_force[i];
                normalize_v2(&force);
                force *= @splat(2, @as(f32, 0.2));
                force += balls_pos[i];

                var center = norm_to_pixel_v2(balls_pos[i]);
                ray.DrawCircleV(center, norm_to_pixel_f(balls_radius[i]), balls_color[i]);
                ray.DrawLineEx(center, norm_to_pixel_v2(force), 3, ray.RED);
            }

            ray.ClearBackground(ray.DARKPURPLE);
        }
        ray.EndDrawing();
    }

    ray.CloseWindow();
}
