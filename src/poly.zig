//! Polynomial operations for Dilithium.
//! Operations on elements of Z_Q[X]/(X^N + 1).

const std = @import("std");

/// Securely clear memory to prevent compiler optimization
/// Uses volatile pointer to ensure memory write is not optimized away
fn secureClearBytes(ptr: []u8) void {
    const volatile_ptr: [*]volatile u8 = @ptrCast(ptr.ptr);
    for (0..ptr.len) |i| {
        volatile_ptr[i] = 0;
    }
}

const simd = @import("simd.zig");
const params = @import("params.zig");
const reduce = @import("reduce.zig");
const rounding = @import("rounding.zig");
const ntt_mod = @import("ntt.zig");
const symmetric = @import("symmetric.zig");

const N = params.N;
const Q = params.Q;
const D = params.D;
const ETA = params.ETA;
const TAU = params.TAU;
const GAMMA1 = params.GAMMA1;
const GAMMA2 = params.GAMMA2;
const SEEDBYTES = params.SEEDBYTES;
const CRHBYTES = params.CRHBYTES;
const CTILDEBYTES = params.CTILDEBYTES;
const POLYZ_PACKEDBYTES = params.POLYZ_PACKEDBYTES;
const POLYETA_PACKEDBYTES = params.POLYETA_PACKEDBYTES;
const POLYT0_PACKEDBYTES = params.POLYT0_PACKEDBYTES;
const POLYT1_PACKEDBYTES = params.POLYT1_PACKEDBYTES;
const POLYW1_PACKEDBYTES = params.POLYW1_PACKEDBYTES;

const STREAM128_BLOCKBYTES = symmetric.STREAM128_BLOCKBYTES;
const STREAM256_BLOCKBYTES = symmetric.STREAM256_BLOCKBYTES;
const SHAKE256_RATE = symmetric.SHAKE256_RATE;

const montgomeryReduce = reduce.montgomeryReduce;
const reduce32 = reduce.reduce32;
const caddq = reduce.caddq;
const power2round = rounding.power2round;
const decompose = rounding.decompose;
const makeHint = rounding.makeHint;
const useHint = rounding.useHint;

/// Polynomial in Z_Q[X]/(X^N + 1).
pub const Poly = struct {
    coeffs: [N]i32 = [_]i32{0} ** N,

    const Self = @This();

    // ==================== Basic Operations ====================

    /// Inplace reduction of all coefficients to representative in [-6283008, 6283008].
    pub fn reduce(self: *Self) void {
        for (&self.coeffs) |*c| {
            c.* = reduce32(c.*);
        }
    }

    /// For all coefficients, add Q if coefficient is negative.
    pub fn caddQ(self: *Self) void {
        for (&self.coeffs) |*c| {
            c.* = caddq(c.*);
        }
    }

    /// Add polynomials: self = a + b. No modular reduction.
    pub fn add(self: *Self, a: *const Self, b: *const Self) void {
        for (&self.coeffs, a.coeffs, b.coeffs) |*c, ca, cb| {
            c.* = ca + cb;
        }
    }

    /// Subtract polynomials: self = a - b. No modular reduction.
    pub fn sub(self: *Self, a: *const Self, b: *const Self) void {
        for (&self.coeffs, a.coeffs, b.coeffs) |*c, ca, cb| {
            c.* = ca - cb;
        }
    }

    /// Multiply polynomial by 2^D without modular reduction.
    /// Assumes input coefficients less than 2^{31-D} in absolute value.
    pub fn shiftl(self: *Self) void {
        for (&self.coeffs) |*c| {
            c.* <<= D;
        }
    }

    // ==================== NTT Operations ====================

    /// Inplace forward NTT. Coefficients can grow by 8*Q in absolute value.
    pub fn ntt(self: *Self) void {
        ntt_mod.ntt(&self.coeffs);
    }

    /// Inplace inverse NTT and multiplication by 2^{32}.
    /// Input coefficients need to be less than Q in absolute value.
    /// Output coefficients are bounded by Q.
    pub fn invnttTomont(self: *Self) void {
        ntt_mod.invnttTomont(&self.coeffs);
    }

    /// Pointwise multiplication in NTT domain with Montgomery reduction.
    /// self = a * b * 2^{-32}
    /// Uses AVX2 acceleration when available.
    pub fn pointwiseMontgomery(self: *Self, a: *const Self, b: *const Self) void {
        if (simd.has_simd) {
            const reduce_avx2 = @import("reduce_avx2.zig");
            var i: usize = 0;
            while (i < N) : (i += 8) {
                const va = simd.load(@ptrCast(&a.coeffs[i]));
                const vb = simd.load(@ptrCast(&b.coeffs[i]));
                const vc = reduce_avx2.pointwiseMont8(va, vb);
                simd.store(@ptrCast(&self.coeffs[i]), vc);
            }
        } else {
            for (&self.coeffs, a.coeffs, b.coeffs) |*c, ca, cb| {
                c.* = montgomeryReduce(@as(i64, ca) * @as(i64, cb));
            }
        }
    }

    // ==================== Rounding Operations ====================

    /// For all coefficients c, compute c0, c1 such that c mod Q = c1*2^D + c0
    /// with -2^{D-1} < c0 <= 2^{D-1}.
    pub fn power2Round(a1: *Self, a0: *Self, a: *const Self) void {
        for (&a1.coeffs, &a0.coeffs, a.coeffs) |*c1, *c0, c| {
            const result = power2round(c);
            c0.* = result.a0;
            c1.* = result.a1;
        }
    }

    /// For all coefficients c, compute high and low bits c0, c1 such that
    /// c mod Q = c1*ALPHA + c0 with -ALPHA/2 < c0 <= ALPHA/2.
    pub fn polyDecompose(a1: *Self, a0: *Self, a: *const Self) void {
        for (&a1.coeffs, &a0.coeffs, a.coeffs) |*c1, *c0, c| {
            const result = decompose(c);
            c0.* = result.a0;
            c1.* = result.a1;
        }
    }

    /// Compute hint polynomial. Returns number of 1 bits.
    pub fn polyMakeHintScalar(h: *Self, a0: *const Self, a1: *const Self) u32 {
        var s: u32 = 0;
        for (&h.coeffs, a0.coeffs, a1.coeffs) |*hc, c0, c1| {
            const hint = makeHint(c0, c1);
            const hint_val: i32 = if (hint) 1 else 0;
            hc.* = hint_val;
            s += @intCast(hint_val);
        }
        return s;
    }
    /// Compute hint polynomial. Returns number of 1 bits.
    /// Compute hint polynomial (SIMD). Returns number of 1 bits.
    pub fn polyMakeHint(h: *Poly, a0: *const Poly, a1: *const Poly) u32 {
        if (!simd.has_simd) {
            return polyMakeHintScalar(h, a0, a1);
        }
        const gamma2_pos: simd.I32x8 = @splat(GAMMA2);
        const gamma2_neg: simd.I32x8 = @splat(-GAMMA2);
        const zero: simd.I32x8 = @splat(0);
        const one: simd.I32x8 = @splat(1);
        const Vec8u = @Vector(8, u32);
        var count: Vec8u = @splat(0);

        const h_ptr: [*]i32 = &h.coeffs;
        const a0_ptr: [*]const i32 = &a0.coeffs;
        const a1_ptr: [*]const i32 = &a1.coeffs;

        comptime var i: usize = 0;
        inline while (i < N) : (i += 8) {
            const c0: simd.I32x8 = a0_ptr[i..][0..8].*;
            const c1: simd.I32x8 = a1_ptr[i..][0..8].*;

            // c0 > GAMMA2
            const gt_pos = c0 > gamma2_pos;

            // c0 < -GAMMA2
            const lt_neg = c0 < gamma2_neg;

            // c0 == -GAMMA2 AND c1 != 0
            const eq_neg = c0 == gamma2_neg;
            const c1_nz = c1 != zero;
            const edge_case = eq_neg & c1_nz;

            // Combine conditions
            const hint_mask = gt_pos | lt_neg | edge_case;

            // Convert bool vector to 0/1 integers
            const hint_vals = @select(i32, hint_mask, one, zero);

            // Store hints
            h_ptr[i..][0..8].* = hint_vals;

            // Accumulate count (reinterpret as unsigned for addition)
            count +%= @bitCast(hint_vals);
        }

        // Horizontal sum
        return @reduce(.Add, count);
    }

    /// Correct high bits according to hint.
    pub fn polyUseHint(b: *Self, a: *const Self, h: *const Self) void {
        for (&b.coeffs, a.coeffs, h.coeffs) |*bc, ac, hc| {
            bc.* = useHint(ac, hc != 0);
        }
    }

    /// Check infinity norm against given bound.
    /// Returns true if norm is strictly smaller than B <= (Q-1)/8.
    pub fn chknormScalar(self: *const Self, bound: i32) bool {
        if (bound > (Q - 1) / 8) return false;

        var failed: i32 = 0;
        for (self.coeffs) |c| {
            // Absolute value without leaking sign
            const mask = c >> 31;
            const abs = c - (mask & (2 * c));

            // Check if abs >= bound
            // This is equivalent to: bound - 1 - abs < 0
            // If abs >= bound, bound - 1 - abs < 0 -> sign bit 1 -> -1 (all 1s)
            // If abs < bound, bound - 1 - abs >= 0 -> sign bit 0 -> 0
            failed |= (bound - 1 - abs) >> 31;
        }

        return failed == 0;
    }

    pub fn chknorm(self: *const Self, bound: i32) bool {
        if (bound > (Q - 1) / 8) return false;
        if (simd.has_simd) {
            return chknorm_avx(self, bound);
        } else {
            return chknormScalar(self, bound);
        }
    }

    pub fn chknorm_avx(self: *const Self, bound: i32) bool {
        if (bound > (Q - 1) / 8) return false;
        const bound_vec: simd.I32x8 = @splat(bound - 1);
        const coeffs_ptr: [*]const i32 = &self.coeffs;

        // 4-way unroll for better ILP
        var f0: simd.I32x8 = @splat(0);
        var f1: simd.I32x8 = @splat(0);
        var f2: simd.I32x8 = @splat(0);
        var f3: simd.I32x8 = @splat(0);

        comptime var i: usize = 0;
        inline while (i < N) : (i += 32) {
            inline for (0..4) |j| {
                const c: simd.I32x8 = coeffs_ptr[i + j * 8 ..][0..8].*;
                const mask = c >> @splat(31);
                const abs = (c ^ mask) -% mask;
                const fail = (bound_vec -% abs) >> @splat(31);

                switch (j) {
                    0 => f0 |= fail,
                    1 => f1 |= fail,
                    2 => f2 |= fail,
                    3 => f3 |= fail,
                    else => unreachable,
                }
            }
        }

        const combined = f0 | f1 | f2 | f3;
        return @reduce(.Or, combined) == 0;
    }

    // ==================== Sampling ====================

    /// Sample polynomial with uniformly random coefficients in [0, Q-1]
    /// using SHAKE128(seed || nonce).
    pub fn uniform(self: *Self, seed: *const [SEEDBYTES]u8, nonce: u16) void {
        const NBLOCKS = (768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES;

        var buf: [NBLOCKS * STREAM128_BLOCKBYTES + 2]u8 = undefined;
        var state = symmetric.Stream128State.init(seed, nonce);
        state.squeezeBlocks(&buf, NBLOCKS);

        var buflen: usize = NBLOCKS * STREAM128_BLOCKBYTES;
        var ctr = rejUniform(self.coeffs[0..], buf[0..buflen]);

        while (ctr < N) {
            const off = buflen % 3;
            for (0..off) |i| {
                buf[i] = buf[buflen - off + i];
            }

            state.squeezeBlocks(buf[off..], 1);
            buflen = STREAM128_BLOCKBYTES + off;
            ctr += rejUniform(self.coeffs[ctr..], buf[0..buflen]);
        }
    }

    /// Sample 4 polynomials with uniformly random coefficients in [0, Q-1]
    /// using parallel SHAKE128(seed || nonces).
    pub fn uniform_4x(p0: *Self, p1: *Self, p2: *Self, p3: *Self, seed: *const [SEEDBYTES]u8, nonces: [4]u16) void {
        const NBLOCKS = (768 + STREAM128_BLOCKBYTES - 1) / STREAM128_BLOCKBYTES;
        const BUF_SIZE = NBLOCKS * STREAM128_BLOCKBYTES + 2;

        var buf0: [BUF_SIZE]u8 = undefined;
        var buf1: [BUF_SIZE]u8 = undefined;
        var buf2: [BUF_SIZE]u8 = undefined;
        var buf3: [BUF_SIZE]u8 = undefined;

        var state = symmetric.Stream128x4State.init(seed, nonces);

        state.squeezeBlocks(buf0[0 .. NBLOCKS * STREAM128_BLOCKBYTES], buf1[0 .. NBLOCKS * STREAM128_BLOCKBYTES], buf2[0 .. NBLOCKS * STREAM128_BLOCKBYTES], buf3[0 .. NBLOCKS * STREAM128_BLOCKBYTES]);

        var buflen0: usize = NBLOCKS * STREAM128_BLOCKBYTES;
        var buflen1: usize = NBLOCKS * STREAM128_BLOCKBYTES;
        var buflen2: usize = NBLOCKS * STREAM128_BLOCKBYTES;
        var buflen3: usize = NBLOCKS * STREAM128_BLOCKBYTES;

        var ctr0 = rejUniform(p0.coeffs[0..], buf0[0..buflen0]);
        var ctr1 = rejUniform(p1.coeffs[0..], buf1[0..buflen1]);
        var ctr2 = rejUniform(p2.coeffs[0..], buf2[0..buflen2]);
        var ctr3 = rejUniform(p3.coeffs[0..], buf3[0..buflen3]);

        while (ctr0 < N or ctr1 < N or ctr2 < N or ctr3 < N) {
            // Preserve leftovers
            const off0 = buflen0 % 3;
            for (0..off0) |i| buf0[i] = buf0[buflen0 - off0 + i];

            const off1 = buflen1 % 3;
            for (0..off1) |i| buf1[i] = buf1[buflen1 - off1 + i];

            const off2 = buflen2 % 3;
            for (0..off2) |i| buf2[i] = buf2[buflen2 - off2 + i];

            const off3 = buflen3 % 3;
            for (0..off3) |i| buf3[i] = buf3[buflen3 - off3 + i];

            // We must squeeze into a temporary buffer if we want to write at offset?
            // Stream128x4State.squeezeBlocks writes to slices.
            // We can pass `buf0[off0..][0..STREAM128_BLOCKBYTES]`

            var tmp0: [STREAM128_BLOCKBYTES]u8 = undefined;
            var tmp1: [STREAM128_BLOCKBYTES]u8 = undefined;
            var tmp2: [STREAM128_BLOCKBYTES]u8 = undefined;
            var tmp3: [STREAM128_BLOCKBYTES]u8 = undefined;

            state.squeezeBlocks(&tmp0, &tmp1, &tmp2, &tmp3);

            @memcpy(buf0[off0..][0..STREAM128_BLOCKBYTES], &tmp0);
            @memcpy(buf1[off1..][0..STREAM128_BLOCKBYTES], &tmp1);
            @memcpy(buf2[off2..][0..STREAM128_BLOCKBYTES], &tmp2);
            @memcpy(buf3[off3..][0..STREAM128_BLOCKBYTES], &tmp3);

            if (ctr0 < N) {
                buflen0 = STREAM128_BLOCKBYTES + off0;
                ctr0 += rejUniform(p0.coeffs[ctr0..], buf0[0..buflen0]);
            }
            if (ctr1 < N) {
                buflen1 = STREAM128_BLOCKBYTES + off1;
                ctr1 += rejUniform(p1.coeffs[ctr1..], buf1[0..buflen1]);
            }
            if (ctr2 < N) {
                buflen2 = STREAM128_BLOCKBYTES + off2;
                ctr2 += rejUniform(p2.coeffs[ctr2..], buf2[0..buflen2]);
            }
            if (ctr3 < N) {
                buflen3 = STREAM128_BLOCKBYTES + off3;
                ctr3 += rejUniform(p3.coeffs[ctr3..], buf3[0..buflen3]);
            }
        }
    }

    /// Sample polynomial with uniformly random coefficients in [-ETA, ETA]
    /// using SHAKE256(seed || nonce).
    pub fn uniformEta(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        const NBLOCKS = if (ETA == 2)
            (136 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES
        else
            (227 + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;

        var buf: [NBLOCKS * STREAM256_BLOCKBYTES]u8 = undefined;
        var state = symmetric.Stream256State.init(seed, nonce);
        state.squeezeBlocks(&buf, NBLOCKS);

        const buflen: usize = NBLOCKS * STREAM256_BLOCKBYTES;
        var ctr = rejEta(self.coeffs[0..], buf[0..buflen]);

        while (ctr < N) {
            state.squeezeBlocks(&buf, 1);
            ctr += rejEta(self.coeffs[ctr..], buf[0..STREAM256_BLOCKBYTES]);
        }
    }

    /// Sample polynomial with uniformly random coefficients in [-(GAMMA1-1), GAMMA1].
    pub fn uniformGamma1(self: *Self, seed: *const [CRHBYTES]u8, nonce: u16) void {
        const NBLOCKS = comptime (POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;

        var buf: [NBLOCKS * STREAM256_BLOCKBYTES]u8 = undefined;
        var state = symmetric.Stream256State.init(seed, nonce);
        state.squeezeBlocks(&buf, NBLOCKS);
        self.zUnpack(buf[0..POLYZ_PACKEDBYTES]);
    }

    /// Sample challenge polynomial with TAU nonzero coefficients in {-1, 1}.
    pub fn challenge(self: *Self, seed: *const [CTILDEBYTES]u8) void {
        var state = std.crypto.hash.sha3.Shake256.init(.{});
        state.update(seed);

        var buf: [SHAKE256_RATE]u8 = undefined;
        state.squeeze(&buf);

        var signs: u64 = 0;
        for (0..8) |i| {
            signs |= @as(u64, buf[i]) << @intCast(8 * i);
        }
        var pos: usize = 8;

        secureClearBytes(std.mem.asBytes(&self.coeffs));

        const start: usize = @as(usize, N) - @as(usize, TAU);
        for (start..N) |i| {
            var b: usize = undefined;

            // Efficient rejection sampling - early exit when valid candidate found
            while (true) {
                if (pos >= SHAKE256_RATE) {
                    state.squeeze(&buf);
                    pos = 0;
                }

                b = buf[pos];
                pos += 1;

                if (b <= i) break;
            }

            self.coeffs[i] = self.coeffs[b];
            self.coeffs[b] = 1 - 2 * @as(i32, @intCast(signs & 1));
            signs >>= 1;
        }
    }
    // ==================== Packing ====================

    /// Bit-pack polynomial with coefficients in [-ETA, ETA].
    pub fn etaPack(self: *const Self, r: *[POLYETA_PACKEDBYTES]u8) void {
        if (ETA == 2) {
            for (0..N / 8) |i| {
                const t = [8]u8{
                    @intCast(ETA - self.coeffs[8 * i + 0]),
                    @intCast(ETA - self.coeffs[8 * i + 1]),
                    @intCast(ETA - self.coeffs[8 * i + 2]),
                    @intCast(ETA - self.coeffs[8 * i + 3]),
                    @intCast(ETA - self.coeffs[8 * i + 4]),
                    @intCast(ETA - self.coeffs[8 * i + 5]),
                    @intCast(ETA - self.coeffs[8 * i + 6]),
                    @intCast(ETA - self.coeffs[8 * i + 7]),
                };
                r[3 * i + 0] = (t[0] >> 0) | (t[1] << 3) | (t[2] << 6);
                r[3 * i + 1] = (t[2] >> 2) | (t[3] << 1) | (t[4] << 4) | (t[5] << 7);
                r[3 * i + 2] = (t[5] >> 1) | (t[6] << 2) | (t[7] << 5);
            }
        } else if (ETA == 4) {
            for (0..N / 2) |i| {
                const t0: u8 = @intCast(ETA - self.coeffs[2 * i + 0]);
                const t1: u8 = @intCast(ETA - self.coeffs[2 * i + 1]);
                r[i] = t0 | (t1 << 4);
            }
        }
    }

    /// Unpack polynomial with coefficients in [-ETA, ETA].
    pub fn etaUnpack(self: *Self, a: *const [POLYETA_PACKEDBYTES]u8) void {
        if (ETA == 2) {
            for (0..N / 8) |i| {
                self.coeffs[8 * i + 0] = @intCast((a[3 * i + 0] >> 0) & 7);
                self.coeffs[8 * i + 1] = @intCast((a[3 * i + 0] >> 3) & 7);
                self.coeffs[8 * i + 2] = @intCast(((a[3 * i + 0] >> 6) | (a[3 * i + 1] << 2)) & 7);
                self.coeffs[8 * i + 3] = @intCast((a[3 * i + 1] >> 1) & 7);
                self.coeffs[8 * i + 4] = @intCast((a[3 * i + 1] >> 4) & 7);
                self.coeffs[8 * i + 5] = @intCast(((a[3 * i + 1] >> 7) | (a[3 * i + 2] << 1)) & 7);
                self.coeffs[8 * i + 6] = @intCast((a[3 * i + 2] >> 2) & 7);
                self.coeffs[8 * i + 7] = @intCast((a[3 * i + 2] >> 5) & 7);

                inline for (0..8) |j| {
                    self.coeffs[8 * i + j] = ETA - self.coeffs[8 * i + j];
                }
            }
        } else if (ETA == 4) {
            for (0..N / 2) |i| {
                self.coeffs[2 * i + 0] = @intCast(a[i] & 0x0F);
                self.coeffs[2 * i + 1] = @intCast(a[i] >> 4);
                self.coeffs[2 * i + 0] = ETA - self.coeffs[2 * i + 0];
                self.coeffs[2 * i + 1] = ETA - self.coeffs[2 * i + 1];
            }
        }
    }

    /// Bit-pack polynomial t1 with coefficients fitting in 10 bits.
    pub fn t1Pack(self: *const Self, r: *[POLYT1_PACKEDBYTES]u8) void {
        for (0..N / 4) |i| {
            r[5 * i + 0] = @truncate(@as(u32, @bitCast(self.coeffs[4 * i + 0])) >> 0);
            r[5 * i + 1] = @truncate((@as(u32, @bitCast(self.coeffs[4 * i + 0])) >> 8) |
                (@as(u32, @bitCast(self.coeffs[4 * i + 1])) << 2));
            r[5 * i + 2] = @truncate((@as(u32, @bitCast(self.coeffs[4 * i + 1])) >> 6) |
                (@as(u32, @bitCast(self.coeffs[4 * i + 2])) << 4));
            r[5 * i + 3] = @truncate((@as(u32, @bitCast(self.coeffs[4 * i + 2])) >> 4) |
                (@as(u32, @bitCast(self.coeffs[4 * i + 3])) << 6));
            r[5 * i + 4] = @truncate(@as(u32, @bitCast(self.coeffs[4 * i + 3])) >> 2);
        }
    }

    /// Unpack polynomial t1 with 10-bit coefficients.
    pub fn t1Unpack(self: *Self, a: *const [POLYT1_PACKEDBYTES]u8) void {
        for (0..N / 4) |i| {
            self.coeffs[4 * i + 0] = @intCast(((@as(u32, a[5 * i + 0]) >> 0) |
                (@as(u32, a[5 * i + 1]) << 8)) & 0x3FF);
            self.coeffs[4 * i + 1] = @intCast(((@as(u32, a[5 * i + 1]) >> 2) |
                (@as(u32, a[5 * i + 2]) << 6)) & 0x3FF);
            self.coeffs[4 * i + 2] = @intCast(((@as(u32, a[5 * i + 2]) >> 4) |
                (@as(u32, a[5 * i + 3]) << 4)) & 0x3FF);
            self.coeffs[4 * i + 3] = @intCast(((@as(u32, a[5 * i + 3]) >> 6) |
                (@as(u32, a[5 * i + 4]) << 2)) & 0x3FF);
        }
    }

    /// Bit-pack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
    pub fn t0Pack(self: *const Self, r: *[POLYT0_PACKEDBYTES]u8) void {
        for (0..N / 8) |i| {
            var t: [8]u32 = undefined;
            inline for (0..8) |j| {
                t[j] = @bitCast((1 << (D - 1)) - self.coeffs[8 * i + j]);
            }

            r[13 * i + 0] = @truncate(t[0]);
            r[13 * i + 1] = @truncate((t[0] >> 8) | (t[1] << 5));
            r[13 * i + 2] = @truncate(t[1] >> 3);
            r[13 * i + 3] = @truncate((t[1] >> 11) | (t[2] << 2));
            r[13 * i + 4] = @truncate((t[2] >> 6) | (t[3] << 7));
            r[13 * i + 5] = @truncate(t[3] >> 1);
            r[13 * i + 6] = @truncate((t[3] >> 9) | (t[4] << 4));
            r[13 * i + 7] = @truncate(t[4] >> 4);
            r[13 * i + 8] = @truncate((t[4] >> 12) | (t[5] << 1));
            r[13 * i + 9] = @truncate((t[5] >> 7) | (t[6] << 6));
            r[13 * i + 10] = @truncate(t[6] >> 2);
            r[13 * i + 11] = @truncate((t[6] >> 10) | (t[7] << 3));
            r[13 * i + 12] = @truncate(t[7] >> 5);
        }
    }

    /// Unpack polynomial t0 with coefficients in ]-2^{D-1}, 2^{D-1}].
    pub fn t0Unpack(self: *Self, a: *const [POLYT0_PACKEDBYTES]u8) void {
        for (0..N / 8) |i| {
            self.coeffs[8 * i + 0] = @intCast((@as(u32, a[13 * i + 0]) |
                (@as(u32, a[13 * i + 1]) << 8)) & 0x1FFF);
            self.coeffs[8 * i + 1] = @intCast(((@as(u32, a[13 * i + 1]) >> 5) |
                (@as(u32, a[13 * i + 2]) << 3) |
                (@as(u32, a[13 * i + 3]) << 11)) & 0x1FFF);
            self.coeffs[8 * i + 2] = @intCast(((@as(u32, a[13 * i + 3]) >> 2) |
                (@as(u32, a[13 * i + 4]) << 6)) & 0x1FFF);
            self.coeffs[8 * i + 3] = @intCast(((@as(u32, a[13 * i + 4]) >> 7) |
                (@as(u32, a[13 * i + 5]) << 1) |
                (@as(u32, a[13 * i + 6]) << 9)) & 0x1FFF);
            self.coeffs[8 * i + 4] = @intCast(((@as(u32, a[13 * i + 6]) >> 4) |
                (@as(u32, a[13 * i + 7]) << 4) |
                (@as(u32, a[13 * i + 8]) << 12)) & 0x1FFF);
            self.coeffs[8 * i + 5] = @intCast(((@as(u32, a[13 * i + 8]) >> 1) |
                (@as(u32, a[13 * i + 9]) << 7)) & 0x1FFF);
            self.coeffs[8 * i + 6] = @intCast(((@as(u32, a[13 * i + 9]) >> 6) |
                (@as(u32, a[13 * i + 10]) << 2) |
                (@as(u32, a[13 * i + 11]) << 10)) & 0x1FFF);
            self.coeffs[8 * i + 7] = @intCast(((@as(u32, a[13 * i + 11]) >> 3) |
                (@as(u32, a[13 * i + 12]) << 5)) & 0x1FFF);

            inline for (0..8) |j| {
                self.coeffs[8 * i + j] = (1 << (D - 1)) - self.coeffs[8 * i + j];
            }
        }
    }

    /// Bit-pack polynomial z with coefficients in [-(GAMMA1-1), GAMMA1].
    pub fn zPack(self: *const Self, r: *[POLYZ_PACKEDBYTES]u8) void {
        if (GAMMA1 == (1 << 17)) {
            for (0..N / 4) |i| {
                var t: [4]u32 = undefined;
                inline for (0..4) |j| {
                    t[j] = @bitCast(GAMMA1 - self.coeffs[4 * i + j]);
                }

                r[9 * i + 0] = @truncate(t[0]);
                r[9 * i + 1] = @truncate(t[0] >> 8);
                r[9 * i + 2] = @truncate((t[0] >> 16) | (t[1] << 2));
                r[9 * i + 3] = @truncate(t[1] >> 6);
                r[9 * i + 4] = @truncate((t[1] >> 14) | (t[2] << 4));
                r[9 * i + 5] = @truncate(t[2] >> 4);
                r[9 * i + 6] = @truncate((t[2] >> 12) | (t[3] << 6));
                r[9 * i + 7] = @truncate(t[3] >> 2);
                r[9 * i + 8] = @truncate(t[3] >> 10);
            }
        } else if (GAMMA1 == (1 << 19)) {
            for (0..N / 2) |i| {
                const t0: u32 = @bitCast(GAMMA1 - self.coeffs[2 * i + 0]);
                const t1: u32 = @bitCast(GAMMA1 - self.coeffs[2 * i + 1]);

                r[5 * i + 0] = @truncate(t0);
                r[5 * i + 1] = @truncate(t0 >> 8);
                r[5 * i + 2] = @truncate((t0 >> 16) | (t1 << 4));
                r[5 * i + 3] = @truncate(t1 >> 4);
                r[5 * i + 4] = @truncate(t1 >> 12);
            }
        }
    }

    /// Unpack polynomial z with coefficients in [-(GAMMA1-1), GAMMA1].
    pub fn zUnpack(self: *Self, a: *const [POLYZ_PACKEDBYTES]u8) void {
        if (GAMMA1 == (1 << 17)) {
            for (0..N / 4) |i| {
                self.coeffs[4 * i + 0] = @intCast((@as(u32, a[9 * i + 0]) |
                    (@as(u32, a[9 * i + 1]) << 8) |
                    (@as(u32, a[9 * i + 2]) << 16)) & 0x3FFFF);
                self.coeffs[4 * i + 1] = @intCast(((@as(u32, a[9 * i + 2]) >> 2) |
                    (@as(u32, a[9 * i + 3]) << 6) |
                    (@as(u32, a[9 * i + 4]) << 14)) & 0x3FFFF);
                self.coeffs[4 * i + 2] = @intCast(((@as(u32, a[9 * i + 4]) >> 4) |
                    (@as(u32, a[9 * i + 5]) << 4) |
                    (@as(u32, a[9 * i + 6]) << 12)) & 0x3FFFF);
                self.coeffs[4 * i + 3] = @intCast(((@as(u32, a[9 * i + 6]) >> 6) |
                    (@as(u32, a[9 * i + 7]) << 2) |
                    (@as(u32, a[9 * i + 8]) << 10)) & 0x3FFFF);

                inline for (0..4) |j| {
                    self.coeffs[4 * i + j] = GAMMA1 - self.coeffs[4 * i + j];
                }
            }
        } else if (GAMMA1 == (1 << 19)) {
            for (0..N / 2) |i| {
                self.coeffs[2 * i + 0] = @intCast((@as(u32, a[5 * i + 0]) |
                    (@as(u32, a[5 * i + 1]) << 8) |
                    (@as(u32, a[5 * i + 2]) << 16)) & 0xFFFFF);
                self.coeffs[2 * i + 1] = @intCast((@as(u32, a[5 * i + 2]) >> 4) |
                    (@as(u32, a[5 * i + 3]) << 4) |
                    (@as(u32, a[5 * i + 4]) << 12));

                self.coeffs[2 * i + 0] = GAMMA1 - self.coeffs[2 * i + 0];
                self.coeffs[2 * i + 1] = GAMMA1 - self.coeffs[2 * i + 1];
            }
        }
    }

    /// Bit-pack polynomial w1 with coefficients in [0,15] or [0,43].
    pub fn w1Pack(self: *const Self, r: *[POLYW1_PACKEDBYTES]u8) void {
        if (GAMMA2 == (Q - 1) / 88) {
            for (0..N / 4) |i| {
                r[3 * i + 0] = @intCast(self.coeffs[4 * i + 0] |
                    (self.coeffs[4 * i + 1] << 6));
                r[3 * i + 1] = @intCast((self.coeffs[4 * i + 1] >> 2) |
                    (self.coeffs[4 * i + 2] << 4));
                r[3 * i + 2] = @intCast((self.coeffs[4 * i + 2] >> 4) |
                    (self.coeffs[4 * i + 3] << 2));
            }
        } else if (GAMMA2 == (Q - 1) / 32) {
            for (0..N / 2) |i| {
                r[i] = @intCast(self.coeffs[2 * i + 0] | (self.coeffs[2 * i + 1] << 4));
            }
        }
    }
};

// ==================== Private Helper Functions ====================

/// Rejection sampling for uniform coefficients in [0, Q-1].
fn rejUniform(a: []i32, buf: []const u8) usize {
    var ctr: usize = 0;
    var pos: usize = 0;

    while (ctr < a.len and pos + 3 <= buf.len) {
        var t: u32 = buf[pos];
        t |= @as(u32, buf[pos + 1]) << 8;
        t |= @as(u32, buf[pos + 2]) << 16;
        t &= 0x7FFFFF;
        pos += 3;

        if (t < Q) {
            a[ctr] = @intCast(t);
            ctr += 1;
        }
    }
    return ctr;
}

/// Rejection sampling for coefficients in [-ETA, ETA].
fn rejEta(a: []i32, buf: []const u8) usize {
    var ctr: usize = 0;
    var pos: usize = 0;

    while (ctr < a.len and pos < buf.len) {
        const t0: u32 = buf[pos] & 0x0F;
        const t1: u32 = buf[pos] >> 4;
        pos += 1;

        if (ETA == 2) {
            if (t0 < 15) {
                const v = t0 -% ((205 * t0) >> 10) * 5;
                a[ctr] = @intCast(2 -% v);
                ctr += 1;
            }
        } else if (ETA == 4) {
            if (t0 < 9) {
                a[ctr] = 4 - @as(i32, @intCast(t0));
                ctr += 1;
            }
        }

        if (ctr < a.len) {
            if (ETA == 2) {
                if (t1 < 15) {
                    const v = t1 -% ((205 * t1) >> 10) * 5;
                    a[ctr] = @intCast(2 -% v);
                    ctr += 1;
                }
            } else if (ETA == 4) {
                if (t1 < 9) {
                    a[ctr] = 4 - @as(i32, @intCast(t1));
                    ctr += 1;
                }
            }
        }
    }
    return ctr;
}

// ==================== Tests ====================

test "poly basic operations" {
    var a = Poly{};
    var b = Poly{};
    var c = Poly{};

    for (&a.coeffs, 0..) |*coeff, i| {
        coeff.* = @intCast(i);
    }
    for (&b.coeffs, 0..) |*coeff, i| {
        coeff.* = @intCast(i * 2);
    }

    c.add(&a, &b);
    try std.testing.expectEqual(@as(i32, 3), c.coeffs[1]);

    c.sub(&b, &a);
    try std.testing.expectEqual(@as(i32, 1), c.coeffs[1]);
}

test "poly chknorm" {
    var a = Poly{};
    for (&a.coeffs) |*c| {
        c.* = 100;
    }
    try std.testing.expect(a.chknorm(200));
    try std.testing.expect(!a.chknorm(50));
}
