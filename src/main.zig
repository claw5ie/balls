const std = @import("std");
const ray = @cImport({
    @cInclude("raylib.h");
    @cInclude("raygui.h");
});

var prng = std.rand.DefaultPrng.init(42);
var rng = prng.random();

const Vec2 = @Vector(2, f32);

const Rect = struct {
    pos: Vec2,
    size: Vec2,
};

const dt: f32 = 1.0 / 60.0;

const SCREEN_WIDTH = 800.0;
const SCREEN_HEIGHT = 600.0;

const UP = SCREEN_HEIGHT / SCREEN_WIDTH;
const DOWN = -SCREEN_HEIGHT / SCREEN_WIDTH;

const Spring = struct {
    start: usize,
    end: usize,

    stiffness: f32,
    rest_length: f32,
    min_length: f32,
    max_length: f32,
    damping_factor: f32,

    half_height: f32,
    hand_length: f32,
    blade_count: usize,
};

const Balls = struct {
    pos: [*]Vec2,
    vel: [*]Vec2,
    force: [*]Vec2,
    disp: [*]Vec2,
    mass: [*]f32,
    radius: [*]f32,
    color: [*]ray.Color,
    edges: []Spring,
    count: usize,
};

fn norm_coord_in_pixels_v2(pos: Vec2) ray.Vector2 {
    return ray.Vector2{
        .x = (pos[0] + 1) * SCREEN_WIDTH / 2,
        .y = (-pos[1] * SCREEN_WIDTH + SCREEN_HEIGHT) / 2,
    };
}

fn norm_coord_in_pixels_r(rect: Rect) ray.Rectangle {
    var pos = norm_coord_in_pixels_v2(rect.pos);

    return ray.Rectangle{
        .x = pos.x,
        .y = pos.y,
        .width = length_in_pixels(rect.size[0]),
        .height = length_in_pixels(rect.size[1]),
    };
}

fn pixels_in_norm_coord(pos: ray.Vector2) Vec2 {
    return Vec2{
        pos.x * 2 / SCREEN_WIDTH - 1,
        (-pos.y * 2 + SCREEN_HEIGHT) / SCREEN_WIDTH,
    };
}

fn length_in_pixels(pos: f32) f32 {
    return pos * SCREEN_WIDTH / 2;
}

fn dist(v0: Vec2, v1: Vec2) f32 {
    var dx = v0[0] - v1[0];
    var dy = v0[1] - v1[1];
    return abs_v2(Vec2{ dx, dy });
}

fn dist_square(v0: Vec2, v1: Vec2) f32 {
    var dx = v0[0] - v1[0];
    var dy = v0[1] - v1[1];
    return abs_square(Vec2{ dx, dy });
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

fn abs_square(v: Vec2) f32 {
    return v[0] * v[0] + v[1] * v[1];
}

fn rand_range(min: f32, max: f32) f32 {
    return @intToFloat(f32, rng.int(u16)) / @as(f32, std.math.maxInt(u16)) * (max - min) + min;
}

fn normalize(v: Vec2) Vec2 {
    var result = v;
    var length = abs_v2(v);
    result /= @splat(2, length);

    return result;
}

fn cross(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[1] - v0[1] * v1[0];
}

fn dot(v0: Vec2, v1: Vec2) f32 {
    return v0[0] * v1[0] + v0[1] * v1[1];
}

fn reflect(vec: Vec2, dir: Vec2) Vec2 {
    var d = normalize(dir);
    var w = Vec2{ cross(vec, d), dot(d, vec) };

    return Vec2{ cross(d, w), dot(d, w) };
}

fn is_inside_rect(point: Vec2, rect: Rect) bool {
    return rect.pos[0] <= point[0] and point[0] <= rect.pos[0] + rect.size[0] and
        rect.pos[1] <= point[1] and point[1] <= rect.pos[1] + rect.size[1];
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
    var disp_length = abs_v2(disp);
    var force = GRAVITY_CONSTANT * balls.mass[i] * balls.mass[j] / (disp_length * disp_length);

    disp /= @splat(2, disp_length);

    balls.force[i] += disp * @splat(2, force);
    balls.force[j] -= disp * @splat(2, force);
}

fn add_spring_force(balls: *Balls, spring: Spring) void {
    var i = spring.start;
    var j = spring.end;
    var disp = balls.pos[j] - balls.pos[i];
    var disp_length = abs_v2(disp);

    disp /= @splat(2, disp_length);

    var damping = dot(disp, balls.vel[j] - balls.vel[i]) * spring.damping_factor;
    var force = spring.stiffness * (disp_length - spring.rest_length) + damping;

    balls.force[i] += @splat(2, force) * disp;
    balls.force[j] -= @splat(2, force) * disp;

    var overlap = disp_length - spring.max_length;
    if (overlap > 0) {
        overlap /= 2;
        balls.disp[i] += @splat(2, overlap) * disp;
        balls.disp[j] -= @splat(2, overlap) * disp;

        var force_disp_comp = dot(balls.force[i], disp);
        if (force_disp_comp < 0) {
            balls.force[i] -= @splat(2, force_disp_comp) * disp;
        }

        force_disp_comp = dot(balls.force[j], disp);
        if (force_disp_comp > 0) {
            balls.force[j] -= @splat(2, force_disp_comp) * disp;
        }
    }

    overlap = disp_length - spring.min_length;
    if (overlap < 0) {
        overlap /= 2;
        balls.disp[i] += @splat(2, overlap) * disp;
        balls.disp[j] -= @splat(2, overlap) * disp;

        var force_disp_comp = dot(balls.force[i], disp);
        if (force_disp_comp > 0) {
            balls.force[i] -= @splat(2, force_disp_comp) * disp;
        }

        force_disp_comp = dot(balls.force[j], disp);
        if (force_disp_comp < 0) {
            balls.force[j] -= @splat(2, force_disp_comp) * disp;
        }
    }
}

var spring_points_buffer: [32]Vec2 = undefined;

fn draw_spring(balls: *Balls, spring: Spring, color: ray.Color) void {
    var start = balls.pos[spring.start];
    var disp = balls.pos[spring.end] - start;
    var full_length = abs_v2(disp);

    disp /= @splat(2, full_length);

    var peak_count = 2 * spring.blade_count - @boolToInt(spring.blade_count % 2 == 0);
    var point_count = peak_count + 2 * 2;
    var points: []Vec2 = spring_points_buffer[0..point_count];

    points[0] = Vec2{ 0, 0 };
    points[1] = Vec2{ spring.hand_length, 0 };

    points[point_count - 2] = Vec2{ full_length - spring.hand_length, 0 };
    points[point_count - 1] = Vec2{ full_length, 0 };

    {
        var base_length = (full_length - 2 * spring.hand_length) / @intToFloat(f32, peak_count);
        var peak = Vec2{
            points[1][0] + base_length / 2,
            (std.math.clamp(full_length, spring.min_length, spring.max_length) - spring.max_length) / (spring.min_length - spring.max_length) * spring.half_height,
        };

        var i: usize = 0;
        while (i < peak_count) : (i += 1) {
            points[i + 2] = peak;
            peak[0] += base_length;
            peak[1] = -peak[1];
        }
    }

    {
        var i: usize = 0;
        while (i < point_count) : (i += 1) {
            var x = points[i][0];
            var y = points[i][1];
            points[i] = .{
                start[0] + disp[0] * x - disp[1] * y,
                start[1] + disp[1] * x + disp[0] * y,
            };
        }
    }

    {
        var i: usize = 1;
        while (i < point_count) : (i += 1) {
            var start_in_pixels = norm_coord_in_pixels_v2(points[i - 1]);
            var end_in_pixels = norm_coord_in_pixels_v2(points[i]);
            ray.DrawLineV(start_in_pixels, end_in_pixels, color);
        }
    }
}

// TODO: make the size of spring hand slighly higher than the radius of ball.
// TODO: add ability to pause simulation, drag balls, and change their properties
//       in runtime.
// TODO: add debugging tools (e.g. toggle physic simulation).

pub fn main() void {
    ray.InitWindow(SCREEN_WIDTH, SCREEN_HEIGHT, "Hello, Traingles!");
    ray.SetTargetFPS(60);

    const MAX_BALL_COUNT = 3;
    var _balls_pos: [MAX_BALL_COUNT]Vec2 = .{
        .{ -0.5, 0.5 },
        .{ -0.5, -0.5 },
        .{ 0.5, -0.5 },
    };
    var _balls_vel: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_force: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_disp: [MAX_BALL_COUNT]Vec2 = undefined;
    var _balls_mass = [MAX_BALL_COUNT]f32{ 2, 5, 10 };
    var _balls_color = [MAX_BALL_COUNT]ray.Color{
        ray.MAGENTA,
        ray.BLUE,
        ray.PINK,
    };
    var _balls_edges = [_]Spring{
        .{
            .start = 0,
            .end = 1,
            .stiffness = 4.0,
            .rest_length = 0.15,
            .min_length = 0.08,
            .max_length = 1.0,
            .damping_factor = 1,
            .half_height = 0.05,
            .hand_length = 0.05,
            .blade_count = 3,
        },
        .{
            .start = 0,
            .end = 2,
            .stiffness = 4.0,
            .rest_length = 0.5,
            .min_length = 0.3,
            .max_length = 0.8,
            .damping_factor = 1,
            .half_height = 0.05,
            .hand_length = 0.05,
            .blade_count = 2,
        },
        .{
            .start = 1,
            .end = 2,
            .stiffness = 4.0,
            .rest_length = 0.8,
            .min_length = 0.3,
            .max_length = 1.3,
            .damping_factor = 1,
            .half_height = 0.05,
            .hand_length = 0.05,
            .blade_count = 1,
        },
    };
    var _balls_radius = [MAX_BALL_COUNT]f32{ 0.05, 0.1, 0.2 };

    var balls = Balls{
        .pos = &_balls_pos,
        .vel = &_balls_vel,
        .force = &_balls_force,
        .disp = &_balls_disp,
        .mass = &_balls_mass,
        .radius = &_balls_radius,
        .color = &_balls_color,
        .edges = _balls_edges[0..],
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

    var ball_menu_rect = Rect{
        .pos = .{ -1, DOWN },
        .size = .{ 0.4, 0.12 },
    };
    var should_draw_ball_menu = false;

    const BallMenuInfo = struct {
        name: [:0]const u8,
        value: *f32,
        upper_limit: f32,
        lower_limit: f32,
    };

    var ball_menu_info = [_]BallMenuInfo{
        .{
            .name = "Mass",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 15,
        },
        .{
            .name = "Radius",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 1.2,
        },
    };

    while (!ray.WindowShouldClose()) {
        ray.BeginDrawing();

        // Clean forces and displacements.
        {
            var i: usize = 0;
            while (i < balls.count) : (i += 1) {
                balls.force[i] = .{ 0, 0 };
                balls.disp[i] = .{ 0, 0 };
            }
        }

        // Apply forces.
        {
            for (balls.edges) |edge| {
                add_spring_force(&balls, edge);
            }
        }

        // Resolve collisions and clip forces.
        {
            {
                var i: usize = 0;
                while (i < balls.count) : (i += 1) {
                    intersect_ball_against_walls(&balls, i);

                    var j: usize = i + 1;
                    while (j < balls.count) : (j += 1) {
                        var disp = balls.pos[j] - balls.pos[i];
                        var disp_length = abs_v2(disp);

                        disp /= @splat(2, disp_length);

                        var overlap = balls.radius[i] + balls.radius[j] - disp_length;
                        if (overlap > 0) {
                            overlap /= 2;

                            balls.disp[i] -= @splat(2, overlap) * disp;
                            balls.disp[j] += @splat(2, overlap) * disp;

                            var force_disp_comp = dot(balls.force[i], disp);
                            if (force_disp_comp > 0) {
                                balls.force[i] -= @splat(2, force_disp_comp) * disp;
                            }

                            if (dot(balls.vel[i], disp) > 0) {
                                balls.vel[i] = reflect(balls.vel[i], Vec2{ -disp[1], disp[0] });
                            }

                            force_disp_comp = dot(balls.force[j], disp);
                            if (force_disp_comp < 0) {
                                balls.force[j] -= @splat(2, force_disp_comp) * disp;
                            }

                            if (dot(balls.vel[j], disp) < 0) {
                                balls.vel[j] = reflect(balls.vel[j], Vec2{ -disp[1], disp[0] });
                            }
                        }
                    }
                }
            }
        }

        // Move balls and draw them.
        {
            {
                var i: usize = 0;
                while (i < balls.count) : (i += 1) {
                    balls.pos[i] += balls.disp[i];

                    balls.vel[i] += balls.force[i] / @splat(2, balls.mass[i]) * @splat(2, dt);
                    balls.pos[i] += balls.vel[i] * @splat(2, dt);
                }
            }

            {
                var i: usize = 0;
                while (i < balls.count) : (i += 1) {
                    var center = norm_coord_in_pixels_v2(balls.pos[i]);
                    ray.DrawCircleV(center, length_in_pixels(balls.radius[i]), balls.color[i]);
                    ray.DrawLineEx(center, norm_coord_in_pixels_v2(balls.pos[i] + balls.force[i] / @splat(2, @as(f32, 4))), 3, ray.RED);
                    ray.DrawLineEx(center, norm_coord_in_pixels_v2(balls.pos[i] + balls.vel[i] / @splat(2, @as(f32, 4))), 3, ray.GREEN);
                }
            }

            for (balls.edges) |edge| {
                draw_spring(&balls, edge, ray.GRAY);
            }
        }

        // Select ball.
        {
            var mouse_pos = pixels_in_norm_coord(ray.GetMousePosition());

            if (ray.IsMouseButtonPressed(ray.MOUSE_BUTTON_LEFT)) {
                if (!is_inside_rect(mouse_pos, ball_menu_rect)) {
                    var is_ball_inside = false;

                    var i: usize = 0;
                    while (i < balls.count) : (i += 1) {
                        var pos = balls.pos[i];
                        pos -= mouse_pos;
                        pos *= pos;

                        var radius = balls.radius[i];
                        radius *= radius;

                        if (pos[0] + pos[1] <= radius) {
                            is_ball_inside = true;
                            ball_menu_info[0].value = &balls.mass[i];
                            ball_menu_info[1].value = &balls.radius[i];
                            break;
                        }
                    }

                    should_draw_ball_menu = is_ball_inside;
                }
            }
        }

        // Draw sliders.
        {
            if (should_draw_ball_menu) {
                var rect = norm_coord_in_pixels_r(ball_menu_rect);
                rect.y -= rect.height;
                rect.height /= @intToFloat(f32, ball_menu_info.len);

                for (ball_menu_info) |slider| {
                    _ = ray.GuiSlider(rect, "", slider.name.ptr, slider.value, slider.lower_limit, slider.upper_limit);
                    rect.y += rect.height;
                }
            }
        }

        ray.ClearBackground(ray.DARKPURPLE);
        ray.EndDrawing();
    }

    ray.CloseWindow();
}
