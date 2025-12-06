//! Polynomial operations for Dilithium.
//! Operations on elements of Z_Q[X]/(X^N + 1).

const std = @import("std");
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
        const simd = @import("simd.zig");
        if (simd.has_avx2) {
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
    pub fn polyMakeHint(h: *Self, a0: *const Self, a1: *const Self) u32 {
        var s: u32 = 0;
        for (&h.coeffs, a0.coeffs, a1.coeffs) |*hc, c0, c1| {
            const hint = makeHint(c0, c1);
            const hint_val: i32 = if (hint) 1 else 0;
            hc.* = hint_val;
            s += @intCast(hint_val);
        }
        return s;
    }
    /// Correct high bits according to hint.
    pub fn polyUseHint(b: *Self, a: *const Self, h: *const Self) void {
        for (&b.coeffs, a.coeffs, h.coeffs) |*bc, ac, hc| {
            bc.* = useHint(ac, hc != 0);
        }
    }

    /// Check infinity norm against given bound.
    /// Returns true if norm is strictly smaller than B <= (Q-1)/8.
    pub fn chknorm(self: *const Self, bound: i32) bool {
        if (bound > (Q - 1) / 8) return false;

        for (self.coeffs) |c| {
            // Absolute value without leaking sign
            const mask = c >> 31;
            const abs = c - (mask & (2 * c));
            if (abs >= bound) return false;
        }
        return true;
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
        const NBLOCKS = (POLYZ_PACKEDBYTES + STREAM256_BLOCKBYTES - 1) / STREAM256_BLOCKBYTES;

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

        @memset(&self.coeffs, 0);

        const start: usize = @as(usize, N) - @as(usize, TAU);
        for (start..N) |i| {
            var b: usize = undefined;
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
            if (t1 < 15 and ctr < a.len) {
                const v = t1 -% ((205 * t1) >> 10) * 5;
                a[ctr] = @intCast(2 -% v);
                ctr += 1;
            }
        } else if (ETA == 4) {
            if (t0 < 9) {
                a[ctr] = 4 - @as(i32, @intCast(t0));
                ctr += 1;
            }
            if (t1 < 9 and ctr < a.len) {
                a[ctr] = 4 - @as(i32, @intCast(t1));
                ctr += 1;
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
