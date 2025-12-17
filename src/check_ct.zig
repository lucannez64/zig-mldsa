const std = @import("std");
const dudect = @import("dudect.zig");
const ntt_mod = @import("ntt.zig");
const poly_mod = @import("poly.zig");
const rounding = @import("rounding.zig");
const params = @import("params.zig");
const sign_mod = @import("sign.zig");
const packing = @import("packing.zig");

const N = params.N;
const GAMMA2 = params.GAMMA2;

fn run_ntt_wrapper(input: *[N]i32) void {
    // We copy input to a buffer to avoid modifying the original during measurement setup if needed,
    // but here we just modify the input in place as ntt expects.
    // However, for dudect, we need to be careful. The measure function measures execution of the function.
    // ntt is in-place.
    var buf = input.*;
    ntt_mod.ntt(&buf);
    // Prevent optimization
    std.mem.doNotOptimizeAway(buf);
}

// Wrapper for chknorm
const ChknormCtx = struct {
    p: poly_mod.Poly,
    bound: i32,
};

fn run_chknorm_wrapper(ctx: *ChknormCtx) void {
    const res = ctx.p.chknorm(ctx.bound);
    std.mem.doNotOptimizeAway(res);
}

// Wrapper for makeHint
const MakeHintCtx = struct {
    a0: i32,
    a1: i32,
};

fn run_make_hint_wrapper(ctx: *MakeHintCtx) void {
    const res = rounding.makeHint(ctx.a0, ctx.a1);
    std.mem.doNotOptimizeAway(res);
}

fn run_sign_wrapper(key: *const packing.SecretKey, msg: []const u8) void {
    const sig = sign_mod.sign(key, msg, "", null) catch unreachable;
    std.mem.doNotOptimizeAway(sig);
}

fn run_keygen_wrapper(seed: *const [params.SEEDBYTES]u8) void {
    const kp = sign_mod.KeyPair.generate(seed);
    std.mem.doNotOptimizeAway(kp);
}

fn run_verify_wrapper(pk: *const packing.PublicKey, msg: []const u8, sig: *const [params.CRYPTO_BYTES]u8) void {
    var valid = true;
    sign_mod.verify(pk, msg, "", sig) catch {
        valid = false;
    };
    std.mem.doNotOptimizeAway(valid);
}

pub fn main() !void {
    // Use crypto random
    const random = std.crypto.random;

    // Test 1: NTT
    {
        var d = dudect.Dudect.init("poly_ntt");

        std.debug.print("Warming up NTT...\n", .{});
        var dummy: [N]i32 = [_]i32{0} ** N;
        for (0..1000) |_| {
            ntt_mod.ntt(&dummy);
        }

        std.debug.print("Starting measurements for {s}...\n", .{d.name});

        const num_measurements = 50000;

        var input0: [N]i32 = [_]i32{0} ** N;
        var input1: [N]i32 = undefined;

        for (0..num_measurements) |i| {
            for (0..N) |j| {
                input1[j] = @as(i32, random.int(i16));
            }

            const class = random.int(u1);

            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_ntt_wrapper, .{&input0});
            } else {
                cycles = dudect.measure_cycles(run_ntt_wrapper, .{&input1});
            }

            d.push(class, cycles);

            if (i % 10000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final NTT t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN NTT!\n", .{});
        }
    }

    // Test 2: Chknorm
    {
        var d = dudect.Dudect.init("poly_chknorm");
        std.debug.print("\nStarting measurements for {s}...\n", .{d.name});

        const bound = (params.Q - 1) / 8; // Max allowed bound

        // Class 0: Fail at index 0
        var ctx0 = ChknormCtx{ .p = .{}, .bound = bound };
        ctx0.p.coeffs[0] = bound + 1;

        // Class 1: Fail at index N-1
        var ctx1 = ChknormCtx{ .p = .{}, .bound = bound };
        ctx1.p.coeffs[N - 1] = bound + 1;

        const num_measurements = 50000;

        for (0..num_measurements) |i| {
            const class = random.int(u1);
            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_chknorm_wrapper, .{&ctx0});
            } else {
                cycles = dudect.measure_cycles(run_chknorm_wrapper, .{&ctx1});
            }
            d.push(class, cycles);
            if (i % 10000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final Chknorm t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN CHKNORM!\n", .{});
        }
    }

    // Test 3: MakeHint
    {
        var d = dudect.Dudect.init("rounding_makeHint");
        std.debug.print("\nStarting measurements for {s}...\n", .{d.name});

        // Class 0: False (a0 = 0)
        var ctx0 = MakeHintCtx{ .a0 = 0, .a1 = 0 };

        // Class 1: True (a0 = GAMMA2 + 1)
        var ctx1 = MakeHintCtx{ .a0 = GAMMA2 + 1, .a1 = 0 };

        const num_measurements = 100000;

        for (0..num_measurements) |i| {
            const class = random.int(u1);
            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_make_hint_wrapper, .{&ctx0});
            } else {
                cycles = dudect.measure_cycles(run_make_hint_wrapper, .{&ctx1});
            }
            d.push(class, cycles);
            if (i % 10000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final MakeHint t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN MAKEHINT!\n", .{});
        }
    }

    // Test 4: Sign
    {
        var d = dudect.Dudect.init("sign_process");
        std.debug.print("\nStarting measurements for {s}...\n", .{d.name});

        // Generate two fixed key pairs for the test
        var seed: [params.SEEDBYTES]u8 = undefined;

        std.crypto.random.bytes(&seed);
        const kp0 = sign_mod.KeyPair.generate(&seed);

        std.crypto.random.bytes(&seed);
        const kp1 = sign_mod.KeyPair.generate(&seed);

        const msg_len = 32;
        var msg: [msg_len]u8 = undefined;

        const num_measurements = 10000;

        for (0..num_measurements) |i| {
            // Random message for each measurement
            std.crypto.random.bytes(&msg);

            const class = random.int(u1);
            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_sign_wrapper, .{ &kp0.secret_key, &msg });
            } else {
                cycles = dudect.measure_cycles(run_sign_wrapper, .{ &kp1.secret_key, &msg });
            }

            d.push(class, cycles);
            if (i % 1000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final Sign t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN SIGN!\n", .{});
        }
    }

    // Test 5: KeyGen
    {
        var d = dudect.Dudect.init("keygen_process");
        std.debug.print("\nStarting measurements for {s}...\n", .{d.name});

        const num_measurements = 10000;
        var seed: [params.SEEDBYTES]u8 = undefined;

        for (0..num_measurements) |i| {
            std.crypto.random.bytes(&seed);
            const class = random.int(u1);

            // Note: KeyGen involves rejection sampling, so it is inherently variable time.
            // We compare two random seeds. The distribution should be identical.
            // If there's a leak based on specific seed bits, it might not show up here
            // unless we target specific classes of seeds.
            // But this checks if "random seed A" vs "random seed B" has different timing characteristics.

            // To make it a proper dudect test, we should have two different input classes.
            // But for KeyGen, any random seed is valid.
            // We'll just run it. The class distinction here is purely nominal (seed A vs seed B).
            // Actually, we should probably generate TWO seeds per iteration, and use one based on class?
            // Yes.

            var seed0: [params.SEEDBYTES]u8 = undefined;
            var seed1: [params.SEEDBYTES]u8 = undefined;
            std.crypto.random.bytes(&seed0);
            std.crypto.random.bytes(&seed1);

            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_keygen_wrapper, .{&seed0});
            } else {
                cycles = dudect.measure_cycles(run_keygen_wrapper, .{&seed1});
            }

            d.push(class, cycles);
            if (i % 1000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final KeyGen t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN KEYGEN (Expected due to rejection sampling)!\n", .{});
        }
    }

    // Test 6: Verify
    {
        var d = dudect.Dudect.init("verify_process");
        std.debug.print("\nStarting measurements for {s}...\n", .{d.name});

        // Generate a fixed key pair and message
        var seed: [params.SEEDBYTES]u8 = undefined;
        std.crypto.random.bytes(&seed);
        const kp = sign_mod.KeyPair.generate(&seed);

        const msg_len = 32;
        var msg: [msg_len]u8 = undefined;
        std.crypto.random.bytes(&msg);

        // Generate a valid signature
        const sig_valid = sign_mod.sign(&kp.secret_key, &msg, "", null) catch unreachable;

        // Generate an invalid signature (flip last byte)
        var sig_invalid = sig_valid;
        sig_invalid[params.CRYPTO_BYTES - 1] ^= 0x01;

        const num_measurements = 10000;

        for (0..num_measurements) |i| {
            // We use fixed key and msg, but vary the signature validity.
            // To avoid branch prediction learning "always valid" or "always invalid", we mix them.
            // However, dudect relies on the fact that we are measuring TWO classes.

            // Ideally we should vary the message too, to ensure we exercise different paths in hashing?
            // But keeping key fixed is fine.
            // Let's vary the message slightly? No, then we need to re-sign.
            // Re-signing is slow.
            // So we stick to fixed message/key, and valid vs invalid signature.

            const class = random.int(u1);
            var cycles: u64 = 0;
            if (class == 0) {
                cycles = dudect.measure_cycles(run_verify_wrapper, .{ &kp.public_key, &msg, &sig_valid });
            } else {
                cycles = dudect.measure_cycles(run_verify_wrapper, .{ &kp.public_key, &msg, &sig_invalid });
            }

            d.push(class, cycles);
            if (i % 1000 == 0 and i > 0) {
                const t = d.t_value();
                std.debug.print("Iter {d}: t-value = {d:.4}\n", .{ i, t });
            }
        }

        const t = d.t_value();
        std.debug.print("Final Verify t-value: {d:.4}\n", .{t});
        if (@abs(t) > 5.0) {
            std.debug.print("POTENTIAL TIMING LEAK IN VERIFY!\n", .{});
        }
    }
}
