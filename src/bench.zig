const std = @import("std");
const builtin = @import("builtin");
const zgg = @import("zgg");

pub fn main() !void {
    std.debug.print("goos: {s}\n", .{@tagName(builtin.os.tag)});
    std.debug.print("goarch: {s}\n", .{@tagName(builtin.cpu.arch)});
    std.debug.print("pkg: zgg\n", .{});

    // Check if SIMD is enabled
    if (zgg.simd.has_simd) {
        std.debug.print("SIMD: enabled\n", .{});
    } else {
        std.debug.print("SIMD: disabled\n", .{});
    }

    const cpu_count = try std.Thread.getCpuCount();

    try benchMlDsa(zgg.MlDsa44, "ML-DSA44", cpu_count);
    try benchMlDsa(zgg.MlDsa65, "ML-DSA65", cpu_count);
    try benchMlDsa(zgg.MlDsa87, "ML-DSA87", cpu_count);
}

fn benchMlDsa(comptime MlDsa: type, name: []const u8, cpu_count: usize) !void {
    // // KeyGen
    // {
    //     var ctx = struct {
    //         pub fn run(_: *const @This()) void {
    //             _ = MlDsa.KeyPair.generate(null);
    //         }
    //     }{};
    //     var buf: [100]u8 = undefined;
    //     const bench_name = try std.fmt.bufPrint(&buf, "BenchmarkKeyGen/{s}", .{name});
    //     try runBenchmark(bench_name, cpu_count, &ctx);
    // }

    // Sign
    {
        var ctx = struct {
            kp: MlDsa.KeyPair,
            msg: []const u8,
            pub fn run(self: *const @This()) void {
                _ = self.kp.secret_key.sign(self.msg, "") catch unreachable;
            }
        }{
            .kp = MlDsa.KeyPair.generate(null),
            .msg = "Hello, World!",
        };
        var buf: [100]u8 = undefined;
        const bench_name = try std.fmt.bufPrint(&buf, "BenchmarkSign/{s}", .{name});
        try runBenchmark(bench_name, cpu_count, &ctx);
    }

    // // Verify
    // {
    //     const kp = MlDsa.KeyPair.generate(null);
    //     const msg = "Hello, World!";
    //     const sig = try kp.secret_key.sign(msg, "");

    //     var ctx = struct {
    //         pk: MlDsa.PublicKey,
    //         msg: []const u8,
    //         sig: @TypeOf(sig),
    //         pub fn run(self: *const @This()) void {
    //             self.pk.verify(self.msg, "", &self.sig) catch unreachable;
    //         }
    //     }{
    //         .pk = kp.public_key,
    //         .msg = msg,
    //         .sig = sig,
    //     };
    //     var buf: [100]u8 = undefined;
    //     const bench_name = try std.fmt.bufPrint(&buf, "BenchmarkVerify/{s}", .{name});
    //     try runBenchmark(bench_name, cpu_count, &ctx);
    // }
}

fn runBenchmark(name: []const u8, cpu_count: usize, ctx: anytype) !void {
    var timer = try std.time.Timer.start();

    var iterations: usize = 1;
    var elapsed: u64 = 0;

    // Aim for 1 second
    const target_ns = 10 * std.time.ns_per_s;

    while (true) {
        timer.reset();
        for (0..iterations) |_| {
            ctx.run();
        }
        elapsed = timer.read();

        if (elapsed >= target_ns) break;

        if (elapsed == 0) {
            iterations *= 10;
        } else {
            const next = iterations * target_ns / elapsed;
            if (next <= iterations) iterations += 1 else iterations = next;
        }
    }

    const ns_per_op = @as(f64, @floatFromInt(elapsed)) / @as(f64, @floatFromInt(iterations));

    // Format: Name-CPUS  Iterations  NsPerOp ns/op  BPerOp B/op  AllocsPerOp allocs/op
    std.debug.print("{s}-{d}\t{d}\t{d:.0} ns/op\t{d} B/op\t{d} allocs/op\n", .{
        name,       cpu_count,
        iterations, ns_per_op,
        0, // B/op
        0, // allocs/op
    });
}
