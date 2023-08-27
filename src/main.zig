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

const Parallelogram = struct {
    pos: Vec2,
    width: Vec2,
    height: Vec2,
};

const dt: f32 = 1.0 / 60.0;

const SCREEN_WIDTH = 800.0;
const SCREEN_HEIGHT = 600.0;

const LEFT = -2.0;
const RIGHT = 2.0;
const DOWN = -SCREEN_HEIGHT / SCREEN_WIDTH * (RIGHT - LEFT) / 2.0;
const UP = -DOWN;

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
    // x in [LEFT, RIGHT] -> x' in [0, SCREEN_WIDTH]
    // y in [DOWN, UP]    -> y' in [SCREEN_HEIGHT, 0]

    return ray.Vector2{
        .x = (pos[0] - LEFT) / (RIGHT - LEFT) * SCREEN_WIDTH,
        .y = (pos[1] - UP) / (DOWN - UP) * SCREEN_HEIGHT,
    };
}

fn norm_coord_in_pixels_r(rect: Rect) ray.Rectangle {
    var pos = norm_coord_in_pixels_v2(rect.pos);
    var result = ray.Rectangle{
        .x = pos.x,
        .y = pos.y,
        .width = length_in_pixels(rect.size[0]),
        .height = length_in_pixels(rect.size[1]),
    };
    result.y -= result.height;
    return result;
}

fn pixels_in_norm_coord(pos: ray.Vector2) Vec2 {
    return Vec2{
        LEFT + pos.x / SCREEN_WIDTH * (RIGHT - LEFT),
        UP + pos.y / SCREEN_HEIGHT * (DOWN - UP),
    };
}

fn length_in_pixels(length: f32) f32 {
    return length * SCREEN_WIDTH / (RIGHT - LEFT);
}

fn dist(v0: Vec2, v1: Vec2) f32 {
    return abs_v2(v0 - v1);
}

fn dist_square(v0: Vec2, v1: Vec2) f32 {
    var result = v0 - v1;
    result *= result;
    return result;
}

fn abs_f32(v: f32) f32 {
    return if (v < 0) -v else v;
}

fn abs_v2(v: Vec2) f32 {
    var u = v;
    var m = @max(abs_f32(v[0]), abs_f32(v[1]));
    if (m == 0) {
        return 0;
    } else {
        u /= @splat(m);
        u *= u;
        return m * std.math.sqrt(u[0] + u[1]);
    }
}

fn rand_range(min: f32, max: f32) f32 {
    return @as(f32, @floatFromInt(rng.int(u16))) / std.math.maxInt(u16) * (max - min) + min;
}

inline fn splat_v2(scalar: f32) Vec2 {
    return @splat(scalar);
}

fn normalize(v: Vec2) Vec2 {
    var result = v;
    var length = abs_v2(v);
    result /= @splat(length);

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

fn clamped_frac(value: f32, lower_bound: f32, upper_bound: f32) f32 {
    return (std.math.clamp(value, lower_bound, upper_bound) - lower_bound) / (upper_bound - lower_bound);
}

fn normal_vector(v: Vec2) Vec2 {
    var n = normalize(v);
    return Vec2{ -n[1], n[0] };
}

fn is_inside_rect(point: Vec2, rect: Rect) bool {
    return rect.pos[0] <= point[0] and point[0] <= rect.pos[0] + rect.size[0] and
        rect.pos[1] <= point[1] and point[1] <= rect.pos[1] + rect.size[1];
}

fn is_inside_parallelogram(point: Vec2, shape: Parallelogram) bool {
    // | width[0] height[0] |^(-1)                        1                           | height[1] -height[0] |
    // | width[1] height[1] |      = ---------------------------------------------- * | -width[1]  width[0]  |
    //                                 width[0] * height[1] - width[1] * height[0]

    var det = shape.width[0] * shape.height[1] - shape.width[1] * shape.height[0];
    std.debug.assert(det > 1e-8);
    var p = (point - shape.pos) / splat_v2(det);

    var x = (shape.height[1] * p[0] - shape.height[0] * p[1]);
    var y = (-shape.width[1] * p[0] + shape.width[0] * p[1]);

    return 0 <= x and x <= 1 and 0 <= y and y <= 1;
}

fn intersect_ball_against_walls(balls: *Balls, index: usize) void {
    var pos = balls.pos[index];
    var radius = balls.radius[index];

    var overlap = pos[0] - radius - LEFT;
    if (overlap < 0) {
        balls.disp[index][0] -= overlap;
        if (balls.force[index][0] < 0) {
            balls.force[index][0] = 0;
        }

        if (balls.vel[index][0] < 0) {
            balls.vel[index][0] = -balls.vel[index][0];
        }
    }

    overlap = pos[0] + radius - RIGHT;
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

    disp /= @splat(disp_length);

    balls.force[i] += disp * splat_v2(force);
    balls.force[j] -= disp * splat_v2(force);
}

fn add_spring_force(balls: *Balls, spring: Spring) void {
    var i = spring.start;
    var j = spring.end;
    var disp = balls.pos[j] - balls.pos[i];
    var disp_length = abs_v2(disp);

    disp /= @splat(disp_length);

    var damping = dot(disp, balls.vel[j] - balls.vel[i]) * spring.damping_factor;
    var force = spring.stiffness * (disp_length - spring.rest_length) + damping;

    balls.force[i] += splat_v2(force) * disp;
    balls.force[j] -= splat_v2(force) * disp;

    var overlap = disp_length - spring.max_length;
    if (overlap > 0) {
        overlap /= 2;
        balls.disp[i] += splat_v2(overlap) * disp;
        balls.disp[j] -= splat_v2(overlap) * disp;

        var force_disp_comp = dot(balls.force[i], disp);
        if (force_disp_comp < 0) {
            balls.force[i] -= splat_v2(force_disp_comp) * disp;
        }

        force_disp_comp = dot(balls.force[j], disp);
        if (force_disp_comp > 0) {
            balls.force[j] -= splat_v2(force_disp_comp) * disp;
        }
    }

    overlap = disp_length - spring.min_length;
    if (overlap < 0) {
        overlap /= 2;
        balls.disp[i] += splat_v2(overlap) * disp;
        balls.disp[j] -= splat_v2(overlap) * disp;

        var force_disp_comp = dot(balls.force[i], disp);
        if (force_disp_comp > 0) {
            balls.force[i] -= splat_v2(force_disp_comp) * disp;
        }

        force_disp_comp = dot(balls.force[j], disp);
        if (force_disp_comp < 0) {
            balls.force[j] -= splat_v2(force_disp_comp) * disp;
        }
    }
}

var spring_points_buffer: [32]Vec2 = undefined;

fn draw_spring(balls: *Balls, spring: Spring, color: ray.Color) void {
    var start = balls.pos[spring.start];
    var disp = balls.pos[spring.end] - start;
    var full_length = abs_v2(disp);

    disp /= @splat(full_length);

    var peak_count = 2 * spring.blade_count - @intFromBool(spring.blade_count % 2 == 0);
    var point_count = peak_count + 2 * 2;
    var points: []Vec2 = spring_points_buffer[0..point_count];

    points[0] = Vec2{ 0, 0 };
    points[1] = Vec2{ spring.hand_length, 0 };

    points[point_count - 2] = Vec2{ full_length - spring.hand_length, 0 };
    points[point_count - 1] = Vec2{ full_length, 0 };

    {
        var base_length = (full_length - 2 * spring.hand_length) / @as(f32, @floatFromInt(peak_count));
        var peak = Vec2{
            points[1][0] + base_length / 2,
            (1.0 - clamped_frac(full_length, spring.min_length, spring.max_length)) * spring.half_height,
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
    var _balls_pos = [_]Vec2{.{ 0, 0 }} ** MAX_BALL_COUNT;
    var _balls_vel = [_]Vec2{.{ 0, 0 }} ** MAX_BALL_COUNT;
    var _balls_force = [_]Vec2{.{ 0, 0 }} ** MAX_BALL_COUNT;
    var _balls_disp = [_]Vec2{.{ 0, 0 }} ** MAX_BALL_COUNT;
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
                rand_range(LEFT + 0.1, RIGHT - 0.1),
                rand_range(DOWN + 0.1, UP - 0.1),
            };
        }
    }

    const SliderInfo = struct {
        name: [:0]const u8,
        value: *f32,
        upper_limit: f32,
        lower_limit: f32,
    };

    var ball_menu_info = [_]SliderInfo{
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

    var spring_menu_info = [_]SliderInfo{
        .{
            .name = "Stiffness",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 10,
        },
        .{
            .name = "Rest Length",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 10,
        },
        .{
            .name = "Minimum Length",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 10,
        },
        .{
            .name = "Maximum Length",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 10,
        },
        .{
            .name = "Damping Factor",
            .value = undefined,
            .lower_limit = 0.01,
            .upper_limit = 10,
        },
    };

    var ball_menu_hitbox = Rect{
        .pos = .{ LEFT, DOWN },
        .size = .{ 1.2, 0.12 * @as(f32, @floatFromInt(ball_menu_info.len)) },
    };
    var should_draw_ball_menu = false;
    var spring_menu_hitbox = Rect{
        .pos = .{ LEFT, DOWN },
        .size = .{ 1.2, 0.12 * @as(f32, @floatFromInt(spring_menu_info.len)) },
    };
    var should_draw_spring_menu = false;

    while (!ray.WindowShouldClose()) {
        ray.BeginDrawing();
        ray.ClearBackground(ray.DARKPURPLE);

        {
            var i: usize = 0;
            while (i < balls.count) : (i += 1) {
                var center = norm_coord_in_pixels_v2(balls.pos[i]);
                ray.DrawCircleV(center, length_in_pixels(balls.radius[i]), balls.color[i]);
                ray.DrawLineEx(center, norm_coord_in_pixels_v2(balls.pos[i] + balls.force[i] / splat_v2(4)), 3, ray.RED);
                ray.DrawLineEx(center, norm_coord_in_pixels_v2(balls.pos[i] + balls.vel[i] / splat_v2(4)), 3, ray.GREEN);
            }
        }

        // Draw springs.
        {
            for (balls.edges) |edge| {
                draw_spring(&balls, edge, ray.GRAY);
            }
        }

        // Draw sliders.
        {
            if (should_draw_ball_menu) {
                var rect = norm_coord_in_pixels_r(ball_menu_hitbox);
                rect.height /= @as(f32, @floatFromInt(ball_menu_info.len));

                for (ball_menu_info) |slider| {
                    _ = ray.GuiSlider(
                        rect,
                        "",
                        slider.name.ptr,
                        slider.value,
                        slider.lower_limit,
                        slider.upper_limit,
                    );
                    rect.y += rect.height;
                }
            }

            if (should_draw_spring_menu) {
                var rect = norm_coord_in_pixels_r(spring_menu_hitbox);
                rect.height /= @as(f32, @floatFromInt(spring_menu_info.len));

                var slider = &spring_menu_info[0];
                _ = ray.GuiSlider(
                    rect,
                    "",
                    slider.name.ptr,
                    slider.value,
                    slider.lower_limit,
                    slider.upper_limit,
                );
                rect.y += rect.height;
                slider = &spring_menu_info[1];
                _ = ray.GuiSlider(
                    rect,
                    "",
                    slider.name.ptr,
                    slider.value,
                    @max(
                        slider.lower_limit,
                        spring_menu_info[2].value.*,
                    ),
                    @min(
                        slider.upper_limit,
                        spring_menu_info[3].value.*,
                    ),
                );
                rect.y += rect.height;
                slider = &spring_menu_info[2];
                _ = ray.GuiSlider(
                    rect,
                    "",
                    slider.name.ptr,
                    slider.value,
                    slider.lower_limit,
                    @min(
                        slider.upper_limit,
                        spring_menu_info[3].value.*,
                    ),
                );
                rect.y += rect.height;
                slider = &spring_menu_info[3];
                _ = ray.GuiSlider(
                    rect,
                    "",
                    slider.name.ptr,
                    slider.value,
                    @max(
                        slider.lower_limit,
                        spring_menu_info[2].value.*,
                    ),
                    slider.upper_limit,
                );
                rect.y += rect.height;
                slider = &spring_menu_info[4];
                _ = ray.GuiSlider(
                    rect,
                    "",
                    slider.name.ptr,
                    slider.value,
                    slider.lower_limit,
                    slider.upper_limit,
                );
            }
        }

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

                        disp /= @splat(disp_length);

                        var overlap = balls.radius[i] + balls.radius[j] - disp_length;
                        if (overlap > 0) {
                            overlap /= 2;

                            balls.disp[i] -= splat_v2(overlap) * disp;
                            balls.disp[j] += splat_v2(overlap) * disp;

                            var force_disp_comp = dot(balls.force[i], disp);
                            if (force_disp_comp > 0) {
                                balls.force[i] -= splat_v2(force_disp_comp) * disp;
                            }

                            if (dot(balls.vel[i], disp) > 0) {
                                balls.vel[i] = reflect(balls.vel[i], Vec2{ -disp[1], disp[0] });
                            }

                            force_disp_comp = dot(balls.force[j], disp);
                            if (force_disp_comp < 0) {
                                balls.force[j] -= splat_v2(force_disp_comp) * disp;
                            }

                            if (dot(balls.vel[j], disp) < 0) {
                                balls.vel[j] = reflect(balls.vel[j], Vec2{ -disp[1], disp[0] });
                            }
                        }
                    }
                }
            }
        }

        // Select ball or spring.
        {
            var mouse_pos = pixels_in_norm_coord(ray.GetMousePosition());

            if (ray.IsMouseButtonPressed(ray.MOUSE_BUTTON_LEFT)) {
                if (!should_draw_spring_menu and !is_inside_rect(mouse_pos, ball_menu_hitbox)) {
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

                if (!should_draw_ball_menu and !is_inside_rect(mouse_pos, spring_menu_hitbox)) {
                    var is_spring_inside = false;

                    for (balls.edges) |*edge| {
                        var width = balls.pos[edge.end] - balls.pos[edge.start];
                        var full_length = abs_v2(width);
                        var height = splat_v2(1.0 - clamped_frac(full_length, edge.min_length, edge.max_length) * edge.half_height * 2) * normal_vector(width);

                        var spring_hitbox = Parallelogram{
                            .pos = balls.pos[edge.start] - height / splat_v2(2),
                            .width = width,
                            .height = height,
                        };

                        if (is_inside_parallelogram(mouse_pos, spring_hitbox)) {
                            is_spring_inside = true;
                            spring_menu_info[0].value = &edge.stiffness;
                            spring_menu_info[1].value = &edge.rest_length;
                            spring_menu_info[2].value = &edge.min_length;
                            spring_menu_info[3].value = &edge.max_length;
                            spring_menu_info[4].value = &edge.damping_factor;
                            break;
                        }
                    }

                    should_draw_spring_menu = is_spring_inside;
                }
            }
        }

        // Move balls
        {
            var i: usize = 0;
            while (i < balls.count) : (i += 1) {
                balls.pos[i] += balls.disp[i];

                balls.vel[i] += balls.force[i] * splat_v2(dt / balls.mass[i]);
                balls.pos[i] += balls.vel[i] * splat_v2(dt);
            }
        }

        ray.EndDrawing();
    }

    ray.CloseWindow();
}
