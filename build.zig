const std = @import("std");

pub fn build(b: *std.Build) void {
    // Option to enable AVX2 optimizations for the library.
    const avx2 = b.option(bool, "avx2", "Enable AVX2 SIMD optimizations") orelse false;

    var target = b.standardTargetOptions(.{});
    if (avx2) {
        target.query.cpu_features_add.addFeature(@intFromEnum(std.Target.x86.Feature.avx2));
    }
    const optimize = b.standardOptimizeOption(.{});

    const mod = b.addModule("zgg", .{
        .root_source_file = b.path("src/root.zig"),
        .target = target,
    });

    const exe = b.addExecutable(.{
        .name = "zgg",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/main.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zgg", .module = mod },
            },
        }),
    });

    b.installArtifact(exe);

    const run_step = b.step("run", "Run the app");
    const run_cmd = b.addRunArtifact(exe);
    run_step.dependOn(&run_cmd.step);

    run_cmd.step.dependOn(b.getInstallStep());

    if (b.args) |args| {
        run_cmd.addArgs(args);
    }

    const mod_tests = b.addTest(.{
        .root_module = mod,
    });

    const run_mod_tests = b.addRunArtifact(mod_tests);

    const exe_tests = b.addTest(.{
        .root_module = exe.root_module,
    });

    const run_exe_tests = b.addRunArtifact(exe_tests);

    const test_step = b.step("test", "Run tests");
    test_step.dependOn(&run_mod_tests.step);
    test_step.dependOn(&run_exe_tests.step);

    // Benchmark step
    const bench_exe = b.addExecutable(.{
        .name = "bench",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/bench.zig"),
            .target = target,
            .optimize = .ReleaseFast, // Benchmarks should be fast
            .imports = &.{
                .{ .name = "zgg", .module = mod },
            },
        }),
    });

    const run_bench = b.addRunArtifact(bench_exe);
    const bench_step = b.step("bench", "Run benchmarks");
    bench_step.dependOn(&run_bench.step);

    const bench_install_step = b.step("bench-install", "Install benchmarks");
    const bench_install_cmd = b.addInstallArtifact(bench_exe, .{});
    bench_install_step.dependOn(&bench_install_cmd.step);

    // Check Constant Time step
    const check_ct_exe = b.addExecutable(.{
        .name = "check_ct",
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/check_ct.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zgg", .module = mod },
            },
        }),
    });

    const run_check_ct = b.addRunArtifact(check_ct_exe);
    const check_ct_step = b.step("check-ct", "Run constant-time checks");
    check_ct_step.dependOn(&run_check_ct.step);

    // NIST KAT step
    const kat_tests = b.addTest(.{
        .root_module = b.createModule(.{
            .root_source_file = b.path("src/nist_kat.zig"),
            .target = target,
            .optimize = optimize,
            .imports = &.{
                .{ .name = "zgg", .module = mod },
            },
        }),
    });

    const run_kat_tests = b.addRunArtifact(kat_tests);
    const check_kat_step = b.step("check-kat", "Run NIST Known Answer Tests");
    check_kat_step.dependOn(&run_kat_tests.step);
}
