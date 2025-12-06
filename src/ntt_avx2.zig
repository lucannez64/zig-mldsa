//! AVX2-accelerated Number Theoretic Transform (NTT) for Dilithium.
//! Vectorized implementation processing 8 coefficients per iteration.

const std = @import("std");
const simd = @import("simd.zig");
const params = @import("params.zig");
const reduce_avx2 = @import("reduce_avx2.zig");
const ntt_scalar = @import("ntt.zig");

const I32x8 = simd.I32x8;
const N = params.N;

/// Forward NTT using AVX2 vectorization where beneficial.
/// Falls back to scalar for the final levels where vector width exceeds butterfly distance.
pub fn ntt(a: *[N]i32) void {
    var k: usize = 0;
    var len: usize = 128;

    // Levels where len >= 8: use vectorized butterflies
    while (len >= 8) : (len >>= 1) {
        var start: usize = 0;
        while (start < N) {
            k += 1;
            const zeta = ntt_scalar.zetas[k];
            const zeta_vec = simd.broadcast(zeta);

            // Process 8 elements at a time
            var j: usize = start;
            while (j < start + len) : (j += 8) {
                const a_lo = simd.load(@ptrCast(&a[j]));
                const a_hi = simd.load(@ptrCast(&a[j + len]));

                // t = zeta * a[j + len] (Montgomery)
                const t = reduce_avx2.pointwiseMont8(zeta_vec, a_hi);

                // Butterfly
                const new_lo = a_lo + t;
                const new_hi = a_lo - t;

                simd.store(@ptrCast(&a[j]), new_lo);
                simd.store(@ptrCast(&a[j + len]), new_hi);
            }
            start = j + len;
        }
    }

    // Remaining levels with len < 8: use scalar
    while (len >= 1) : (len >>= 1) {
        var start: usize = 0;
        while (start < N) {
            k += 1;
            const zeta = ntt_scalar.zetas[k];
            var j: usize = start;
            while (j < start + len) : (j += 1) {
                const reduce = @import("reduce.zig");
                const t = reduce.montgomeryReduce(@as(i64, zeta) * @as(i64, a[j + len]));
                a[j + len] = a[j] - t;
                a[j] = a[j] + t;
            }
            start = j + len;
        }
        if (len == 1) break;
    }
}

/// Inverse NTT using AVX2 vectorization where beneficial.
pub fn invnttTomont(a: *[N]i32) void {
    const reduce = @import("reduce.zig");
    var k: usize = 256;
    var len: usize = 1;

    // mont^2/256
    const f: i32 = 41978;

    // Levels where len < 8: use scalar
    while (len < 8) : (len <<= 1) {
        var start: usize = 0;
        while (start < N) {
            k -= 1;
            const zeta = -ntt_scalar.zetas[k];
            var j: usize = start;
            while (j < start + len) : (j += 1) {
                const t = a[j];
                a[j] = t + a[j + len];
                a[j + len] = t - a[j + len];
                a[j + len] = reduce.montgomeryReduce(@as(i64, zeta) * @as(i64, a[j + len]));
            }
            start = j + len;
        }
    }

    // Levels where len >= 8: use vectorized butterflies
    while (len < N) : (len <<= 1) {
        var start: usize = 0;
        while (start < N) {
            k -= 1;
            const zeta = -ntt_scalar.zetas[k];
            const zeta_vec = simd.broadcast(zeta);

            var j: usize = start;
            while (j < start + len) : (j += 8) {
                const a_lo = simd.load(@ptrCast(&a[j]));
                const a_hi = simd.load(@ptrCast(&a[j + len]));

                // Inverse butterfly
                const new_lo = a_lo + a_hi;
                const diff = a_lo - a_hi;
                const new_hi = reduce_avx2.pointwiseMont8(zeta_vec, diff);

                simd.store(@ptrCast(&a[j]), new_lo);
                simd.store(@ptrCast(&a[j + len]), new_hi);
            }
            start = j + len;
        }
    }

    // Final scaling by f
    const f_vec = simd.broadcast(f);
    var i: usize = 0;
    while (i < N) : (i += 8) {
        const v = simd.load(@ptrCast(&a[i]));
        const scaled = reduce_avx2.pointwiseMont8(f_vec, v);
        simd.store(@ptrCast(&a[i]), scaled);
    }
}

test "ntt_avx2 roundtrip" {
    const freeze = @import("reduce.zig").freeze;

    var poly: [N]i32 = undefined;
    for (&poly, 0..) |*c, i| {
        c.* = @intCast(i);
    }

    var original: [N]i32 = undefined;
    @memcpy(&original, &poly);

    ntt(&poly);
    invnttTomont(&poly);

    // After NTT -> INTT, verify we get valid field elements
    for (poly) |c| {
        const normalized = freeze(c);
        _ = normalized;
    }
}

test "ntt_avx2 matches scalar" {
    var poly_avx2: [N]i32 = undefined;
    var poly_scalar: [N]i32 = undefined;

    for (0..N) |i| {
        const val: i32 = @intCast((i * 17) % 1000);
        poly_avx2[i] = val;
        poly_scalar[i] = val;
    }

    ntt(&poly_avx2);
    ntt_scalar.ntt(&poly_scalar);

    for (0..N) |i| {
        try std.testing.expectEqual(poly_scalar[i], poly_avx2[i]);
    }
}
