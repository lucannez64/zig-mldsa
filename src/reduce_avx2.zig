//! AVX2-accelerated Montgomery reduction and modular arithmetic.
//! Vectorized versions processing 8 elements at a time.

const std = @import("std");
const simd = @import("simd.zig");
const params = @import("params.zig");

const I32x8 = simd.I32x8;
const I64x8 = simd.I64x8;
const Q = params.Q;

// q^(-1) mod 2^32
const QINV: i32 = 58728449;

const Q_VEC: I32x8 = simd.broadcast(Q);
const QINV_VEC: I32x8 = simd.broadcast(QINV);

/// Vectorized Montgomery reduction for 8 elements.
/// For each element a with -2^31 * Q <= a <= Q * 2^31,
/// computes r ≡ a * 2^-32 (mod Q) such that -Q < r < Q.
pub inline fn montgomeryReduce8(a: I64x8) I32x8 {
    // t = (a_lo * QINV) mod 2^32
    const a_lo = simd.truncate64(a);
    const t = simd.mulLo(a_lo, QINV_VEC);

    // product = t * Q (as i64)
    const product = simd.mulWide(t, Q_VEC);

    // result = (a - product) >> 32
    const diff = a - product;
    return simd.shr32(diff);
}

/// Vectorized modular reduction for 8 elements.
/// For elements a with a <= 2^31 - 2^22 - 1,
/// computes r ≡ a (mod Q) such that -6283008 <= r <= 6283008.
pub inline fn reduce32x8(a: I32x8) I32x8 {
    const shift_val: I32x8 = simd.broadcast(1 << 22);
    const t = (a + shift_val) >> @splat(23);
    return a - simd.mulLo(t, Q_VEC);
}

/// Vectorized conditional add Q for 8 elements.
/// Adds Q if input coefficient is negative.
pub inline fn caddq8(a: I32x8) I32x8 {
    const mask = a >> @splat(31);
    return a + (mask & Q_VEC);
}

/// Pointwise Montgomery multiplication of 8 elements.
/// result = a * b * 2^{-32} mod Q
pub inline fn pointwiseMont8(a: I32x8, b: I32x8) I32x8 {
    const product = simd.mulWide(a, b);
    return montgomeryReduce8(product);
}

/// Lazy pointwise multiplication.
/// Performs a * b (low 32 bits) without reduction.
/// Assumes result fits in i32 (which is true for Kyber coefficients * zeta).
pub inline fn pointwiseLazyMul(a: I32x8, b: I32x8) I32x8 {
    return simd.mulLo(a, b);
}

/// Lazy reduction.
/// Reduces only when coefficients exceed a threshold to prevent overflow.
/// Uses Montgomery reduction when reducing.
pub inline fn lazyReduce(x: I32x8) I32x8 {
    // Threshold to prevent overflow in subsequent additions.
    // Lower threshold to 100M to be safe.
    const threshold: I32x8 = @splat(100_000_000);
    const neg_threshold: I32x8 = @splat(-100_000_000);

    // Check if any lane exceeds bounds
    const mask = (x > threshold) | (x < neg_threshold);

    // If mask is all zeros, we can return x (fast path).
    // But for vector selection, we compute reduced version.
    const reduced = montgomeryReduce8(simd.widen(x));

    return @select(i32, mask, reduced, x);
}

test "reduce_avx2 montgomery" {
    const reduce_scalar = @import("reduce.zig");

    // Test that vectorized matches scalar
    const test_vals: [8]i64 = .{ 1000, -2000, 3000, -4000, 5000, -6000, 7000, -8000 };
    const vec_result = montgomeryReduce8(test_vals);
    const result_arr: [8]i32 = vec_result;

    for (0..8) |i| {
        const scalar_result = reduce_scalar.montgomeryReduce(test_vals[i]);
        try std.testing.expectEqual(scalar_result, result_arr[i]);
    }
}

test "reduce_avx2 pointwise" {
    const reduce_scalar = @import("reduce.zig");

    const a: [8]i32 = .{ 100, 200, 300, 400, 500, 600, 700, 800 };
    const b: [8]i32 = .{ 10, 20, 30, 40, 50, 60, 70, 80 };

    const vec_result = pointwiseMont8(a, b);
    const result_arr: [8]i32 = vec_result;

    for (0..8) |i| {
        const scalar_result = reduce_scalar.montgomeryReduce(@as(i64, a[i]) * @as(i64, b[i]));
        try std.testing.expectEqual(scalar_result, result_arr[i]);
    }
}
