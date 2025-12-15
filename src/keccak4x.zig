const std = @import("std");
const simd = @import("simd.zig");
const builtin = @import("builtin");

pub const U64x4 = @Vector(4, u64);

/// Keccak-f[1600] state containing 4 parallel states.
/// 25 lanes, each lane is a 4-element vector of u64.
pub const KeccakF1600x4 = struct {
    s: [25]U64x4 align(32) = [_]U64x4{@splat(0)} ** 25,
    pub fn permute(self: *@This()) void {
        @setEvalBranchQuota(10000);
        var s = &self.s;
        inline for (round_constants) |rc| {
            // Theta
            var c: [5]U64x4 = undefined;
            var d: [5]U64x4 = undefined;

            inline for (0..5) |i| {
                c[i] = s[i] ^ s[i + 5] ^ s[i + 10] ^ s[i + 15] ^ s[i + 20];
            }

            inline for (0..5) |i| {
                d[i] = c[(i + 4) % 5] ^ rol(c[(i + 1) % 5], 1);
            }

            inline for (0..25) |i| {
                s[i] ^= d[i % 5];
            }

            // Rho + Pi
            // We need a temporary buffer because Pi permutes positions
            var b: [25]U64x4 = undefined;
            inline for (0..5) |y| {
                inline for (0..5) |x| {
                    const index = x + 5 * y;
                    b[y + 5 * ((2 * x + 3 * y) % 5)] = rol(s[index], rho_offsets[index]);
                }
            }

            // Chi
            inline for (0..5) |y| {
                const row = y * 5;

                // Load entire row into registers (5 vector loads)
                const b0 = b[row + 0];
                const b1 = b[row + 1];
                const b2 = b[row + 2];
                const b3 = b[row + 3];
                const b4 = b[row + 4];

                // Precompute negations (vectorized across 4 lanes)
                const not_b0 = ~b0;
                const not_b1 = ~b1;
                const not_b2 = ~b2;
                const not_b3 = ~b3;
                const not_b4 = ~b4;

                // Compute Chi with explicit unrolling - no modulo ops
                // Each: 1 vector AND, 1 vector XOR
                s[row + 0] = b0 ^ (not_b1 & b2);
                s[row + 1] = b1 ^ (not_b2 & b3);
                s[row + 2] = b2 ^ (not_b3 & b4);
                s[row + 3] = b3 ^ (not_b4 & b0);
                s[row + 4] = b4 ^ (not_b0 & b1);
            }

            // Iota
            s[0] ^= @as(U64x4, @splat(rc));
        }
    }
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
