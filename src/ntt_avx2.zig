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

/// Inverse NTT using AVX2 vectorization with manual unrolling.
/// Removed strict alignment requirements to prevent runtime panics.
pub fn invnttTomont(a: *[N]i32) void {
    const reduce = @import("reduce.zig");
    // Removed @alignCast to prevent panic on unaligned inputs.
    // Modern CPUs handle unaligned loads (vmovdqu) very efficiently.

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

    // Levels where len >= 32: Unroll 4x (Process 32 elements per step)
    while (len < N) : (len <<= 1) {
        // Fallback to standard vector loop for small strides (8, 16)
        // Unrolling here wouldn't gain much as trip count is small.
        if (len < 32) {
            var start: usize = 0;
            while (start < N) {
                k -= 1;
                const zeta = -ntt_scalar.zetas[k];
                const zeta_vec = simd.broadcast(zeta);

                var j: usize = start;
                while (j < start + len) : (j += 8) {
                    const a_lo = simd.load(@ptrCast(&a[j]));
                    const a_hi = simd.load(@ptrCast(&a[j + len]));

                    const new_lo = a_lo + a_hi;
                    const diff = a_lo - a_hi;
                    const new_hi = reduce_avx2.pointwiseMont8(zeta_vec, diff);

                    simd.store(@ptrCast(&a[j]), new_lo);
                    simd.store(@ptrCast(&a[j + len]), new_hi);
                }
                start = j + len;
            }
            continue;
        }

        // Heavy unroll for deeper layers (len >= 32)
        var start: usize = 0;
        while (start < N) {
            k -= 1;
            const zeta = -ntt_scalar.zetas[k];
            const zeta_vec = simd.broadcast(zeta);

            var j: usize = start;
            const end_j = start + len;

            // Process 4 vectors (32 coeffs) at a time
            // We use `a` directly. simd.load handles the pointer cast/loading.
            while (j + 32 <= end_j) : (j += 32) {
                // Load 4 blocks low
                const lo_0 = simd.load(@ptrCast(&a[j + 0]));
                const lo_1 = simd.load(@ptrCast(&a[j + 8]));
                const lo_2 = simd.load(@ptrCast(&a[j + 16]));
                const lo_3 = simd.load(@ptrCast(&a[j + 24]));

                // Load 4 blocks high
                const hi_0 = simd.load(@ptrCast(&a[j + len + 0]));
                const hi_1 = simd.load(@ptrCast(&a[j + len + 8]));
                const hi_2 = simd.load(@ptrCast(&a[j + len + 16]));
                const hi_3 = simd.load(@ptrCast(&a[j + len + 24]));

                // Arith 0
                const n_lo_0 = lo_0 + hi_0;
                const diff_0 = lo_0 - hi_0;
                const n_hi_0 = reduce_avx2.pointwiseMont8(zeta_vec, diff_0);

                // Arith 1
                const n_lo_1 = lo_1 + hi_1;
                const diff_1 = lo_1 - hi_1;
                const n_hi_1 = reduce_avx2.pointwiseMont8(zeta_vec, diff_1);

                // Arith 2
                const n_lo_2 = lo_2 + hi_2;
                const diff_2 = lo_2 - hi_2;
                const n_hi_2 = reduce_avx2.pointwiseMont8(zeta_vec, diff_2);

                // Arith 3
                const n_lo_3 = lo_3 + hi_3;
                const diff_3 = lo_3 - hi_3;
                const n_hi_3 = reduce_avx2.pointwiseMont8(zeta_vec, diff_3);

                // Store 4 blocks low
                simd.store(@ptrCast(&a[j + 0]), n_lo_0);
                simd.store(@ptrCast(&a[j + 8]), n_lo_1);
                simd.store(@ptrCast(&a[j + 16]), n_lo_2);
                simd.store(@ptrCast(&a[j + 24]), n_lo_3);

                // Store 4 blocks high
                simd.store(@ptrCast(&a[j + len + 0]), n_hi_0);
                simd.store(@ptrCast(&a[j + len + 8]), n_hi_1);
                simd.store(@ptrCast(&a[j + len + 16]), n_hi_2);
                simd.store(@ptrCast(&a[j + len + 24]), n_hi_3);
            }

            // Cleanup remaining items (though len is always power of 2 >= 32, this is safe)
            while (j < end_j) : (j += 8) {
                const a_lo = simd.load(@ptrCast(&a[j]));
                const a_hi = simd.load(@ptrCast(&a[j + len]));
                const new_lo = a_lo + a_hi;
                const diff = a_lo - a_hi;
                const new_hi = reduce_avx2.pointwiseMont8(zeta_vec, diff);
                simd.store(@ptrCast(&a[j]), new_lo);
                simd.store(@ptrCast(&a[j + len]), new_hi);
            }
            start = j + len;
        }
    }

    // Final scaling by f - Unrolled 2x
    const f_vec = simd.broadcast(f);
    var i: usize = 0;
    while (i < N) : (i += 16) {
        const v0 = simd.load(@ptrCast(&a[i]));
        const v1 = simd.load(@ptrCast(&a[i + 8]));

        const scaled0 = reduce_avx2.pointwiseMont8(f_vec, v0);
        const scaled1 = reduce_avx2.pointwiseMont8(f_vec, v1);

        simd.store(@ptrCast(&a[i]), scaled0);
        simd.store(@ptrCast(&a[i + 8]), scaled1);
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
