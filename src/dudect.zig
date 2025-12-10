const std = @import("std");
const builtin = @import("builtin");

pub const Dudect = struct {
    name: []const u8,

    // Statistics for class 1
    n1: f64 = 0,
    mean1: f64 = 0,
    m2_1: f64 = 0, // Sum of squares of differences from the current mean

    // Statistics for class 2
    n2: f64 = 0,
    mean2: f64 = 0,
    m2_2: f64 = 0,

    pub fn init(name: []const u8) Dudect {
        return .{
            .name = name,
        };
    }

    // Welford's online algorithm
    fn update(n: *f64, mean: *f64, m2: *f64, value: f64) void {
        n.* += 1;
        const delta = value - mean.*;
        mean.* += delta / n.*;
        const delta2 = value - mean.*;
        m2.* += delta * delta2;
    }

    pub fn push(self: *Dudect, class: u1, cycles: u64) void {
        const val = @as(f64, @floatFromInt(cycles));
        if (class == 0) {
            update(&self.n1, &self.mean1, &self.m2_1, val);
        } else {
            update(&self.n2, &self.mean2, &self.m2_2, val);
        }
    }

    pub fn t_value(self: Dudect) f64 {
        if (self.n1 < 2 or self.n2 < 2) return 0;

        const var1 = self.m2_1 / (self.n1 - 1);
        const var2 = self.m2_2 / (self.n2 - 1);

        const num = self.mean1 - self.mean2;
        const den = std.math.sqrt((var1 / self.n1) + (var2 / self.n2));

        if (den == 0) return 0;
        return num / den;
    }

    // Returns max t-value, generally > 4.5 is considered a failure (99.999% confidence)
    // but usually we want to see if it grows with N.
    // For automated testing, a threshold like 5.0 is often used.
};

pub fn measure_cycles(func: anytype, args: anytype) u64 {
    const start = rdtsc();
    @call(.auto, func, args);
    const end = rdtsc();
    // Handle overflow (though unlikely to matter for short durations)
    if (end < start) return 0;
    return end - start;
}

fn rdtsc() u64 {
    if (builtin.cpu.arch == .x86 or builtin.cpu.arch == .x86_64) {
        var low: u32 = undefined;
        var high: u32 = undefined;
        asm volatile ("rdtsc"
            : [low] "={eax}" (low),
              [high] "={edx}" (high),
        );
        return (@as(u64, high) << 32) | low;
    } else {
        // Fallback for other architectures (e.g. ARM)
        // Ideally use performace counters, but std.time is a fallback
        return @as(u64, @intCast(std.time.nanoTimestamp()));
    }
}

// Prepare inputs for the test
// The user needs to provide a way to generate inputs for two classes.
// Class 0: Fix vs Random?
// Class 1: Random?
