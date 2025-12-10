const std = @import("std");
const simd = @import("simd.zig");
const builtin = @import("builtin");

pub const U64x4 = @Vector(4, u64);

/// Keccak-f[1600] state containing 4 parallel states.
/// 25 lanes, each lane is a 4-element vector of u64.
pub const State4x = struct {
    s: [25]U64x4,
};

/// Round constants for Keccak-f[1600]
const round_constants = [_]u64{
    0x0000000000000001, 0x0000000000008082, 0x800000000000808a, 0x8000000080008000,
    0x000000000000808b, 0x0000000080000001, 0x8000000080008081, 0x8000000000008009,
    0x000000000000008a, 0x0000000000000088, 0x0000000080008009, 0x000000008000000a,
    0x000000008000808b, 0x800000000000008b, 0x8000000000008089, 0x8000000000008003,
    0x8000000000008002, 0x8000000000000080, 0x000000000000800a, 0x800000008000000a,
    0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008,
};

/// Rotation offsets for Rho step
const rho_offsets = [_]u6{
    0,  1,  62, 28, 27,
    36, 44, 6,  55, 20,
    3,  10, 43, 25, 39,
    41, 45, 15, 21, 8,
    18, 2,  61, 56, 14,
};

inline fn rol(x: U64x4, comptime n: u6) U64x4 {
    if (n == 0) return x;
    if (builtin.cpu.arch == .x86_64 and std.Target.x86.featureSetHas(builtin.cpu.features, .bmi2)) {
        // RORX is faster on some microarchitectures
        // return @bitCast(std.simd.rotateElementsLeft(x, n));
        // std.simd.rotateElementsLeft not available?
    }
    const shift: u64 = n;
    return (x << @splat(shift)) | (x >> @splat(64 - shift));
}

pub fn keccakF1600_4x(state: *State4x) void {
    @setEvalBranchQuota(10000);
    inline for (round_constants) |rc| {
        // Theta
        var c: [5]U64x4 = undefined;
        var d: [5]U64x4 = undefined;

        inline for (0..5) |i| {
            c[i] = state.s[i] ^ state.s[i + 5] ^ state.s[i + 10] ^ state.s[i + 15] ^ state.s[i + 20];
        }

        inline for (0..5) |i| {
            d[i] = c[(i + 4) % 5] ^ rol(c[(i + 1) % 5], 1);
        }

        inline for (0..25) |i| {
            state.s[i] ^= d[i % 5];
        }

        // Rho + Pi
        // We need a temporary buffer because Pi permutes positions
        var b: [25]U64x4 = undefined;
        inline for (0..5) |y| {
            inline for (0..5) |x| {
                const index = x + 5 * y;
                b[y + 5 * ((2 * x + 3 * y) % 5)] = rol(state.s[index], rho_offsets[index]);
            }
        }

        // Chi
        inline for (0..5) |y| {
            inline for (0..5) |x| {
                const index = x + 5 * y;
                state.s[index] = b[index] ^ ((~b[(x + 1) % 5 + 5 * y]) & b[(x + 2) % 5 + 5 * y]);
            }
        }

        // Iota
        state.s[0] ^= @as(U64x4, @splat(rc));
    }
}
