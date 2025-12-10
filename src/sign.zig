//! Dilithium signature scheme.
//! ML-DSA (FIPS 204) implementation.

const std = @import("std");
const assert = std.debug.assert;

/// Securely clear memory to prevent compiler optimization
/// Uses volatile pointer to ensure memory write is not optimized away
fn secureClear(comptime T: type, ptr: *T) void {
    const volatile_ptr: *volatile T = @ptrCast(ptr);
    volatile_ptr.* = std.mem.zeroes(T);
}

/// Securely clear a byte array to prevent compiler optimization
fn secureClearBytes(ptr: []u8) void {
    const volatile_ptr: [*]volatile u8 = @ptrCast(ptr.ptr);
    for (0..ptr.len) |i| {
        volatile_ptr[i] = 0;
    }
}
const params = @import("params.zig");
const poly_mod = @import("poly.zig");
const polyvec_mod = @import("polyvec.zig");
const packing = @import("packing.zig");
const config = @import("config.zig");
const simd = @import("simd.zig");

const K = params.K;
const L = params.L;
const N = params.N;
const OMEGA = params.OMEGA;
const BETA = params.BETA;
const GAMMA1 = params.GAMMA1;
const GAMMA2 = params.GAMMA2;
const SEEDBYTES = params.SEEDBYTES;
const CRHBYTES = params.CRHBYTES;
const TRBYTES = params.TRBYTES;
const RNDBYTES = params.RNDBYTES;
const CTILDEBYTES = params.CTILDEBYTES;
const POLYW1_PACKEDBYTES = params.POLYW1_PACKEDBYTES;
const CRYPTO_PUBLICKEYBYTES = params.CRYPTO_PUBLICKEYBYTES;
const CRYPTO_SECRETKEYBYTES = params.CRYPTO_SECRETKEYBYTES;
const CRYPTO_BYTES = params.CRYPTO_BYTES;

const Poly = poly_mod.Poly;
const PolyVecL = polyvec_mod.PolyVecL;
const PolyVecK = polyvec_mod.PolyVecK;
const Mat = polyvec_mod.Mat;

const Shake256 = std.crypto.hash.sha3.Shake256;

pub const SigningError = error{
    ContextTooLong,
    MessageTooLong,
    InvalidSignatureLength,
    InvalidSignature,
    VerificationFailed,
};

/// Configuration for signing.
pub const SigningMode = enum {
    randomized,
    deterministic,
};

const signing_mode: SigningMode = config.signing_mode;

/// Default signing mode (can be changed at comptime).
/// Dilithium key pair.
pub const KeyPair = struct {
    public_key: packing.PublicKey,
    secret_key: packing.SecretKey,

    const Self = @This();

    /// Generate a new key pair.
    pub fn generate(random: ?*const [SEEDBYTES]u8) Self {
        var seedbuf: [2 * SEEDBYTES + CRHBYTES]u8 = undefined;
        defer secureClearBytes(&seedbuf); // Zero sensitive seed material

        // Get randomness for rho, rhoprime and key
        if (random) |r| {
            @memcpy(seedbuf[0..SEEDBYTES], r);
        } else {
            std.crypto.random.bytes(seedbuf[0..SEEDBYTES]);
        }

        seedbuf[SEEDBYTES + 0] = K;
        seedbuf[SEEDBYTES + 1] = L;

        var hasher = Shake256.init(.{});
        hasher.update(seedbuf[0 .. SEEDBYTES + 2]);
        hasher.squeeze(seedbuf[0 .. 2 * SEEDBYTES + CRHBYTES]);

        const rho = seedbuf[0..SEEDBYTES];
        const rhoprime = seedbuf[SEEDBYTES..][0..CRHBYTES];
        const key = seedbuf[SEEDBYTES + CRHBYTES ..][0..SEEDBYTES];

        // Expand matrix
        var mat: Mat = .{};
        mat.expand(rho);

        // Sample short vectors s1 and s2
        var s1: PolyVecL = .{};
        var s2: PolyVecK = .{};
        s1.uniformEta(rhoprime, 0);
        s2.uniformEta(rhoprime, L);

        // Matrix-vector multiplication
        var s1hat = s1;
        s1hat.ntt();

        var t1: PolyVecK = .{};
        Mat.pointwiseMontgomery(&t1, &mat, &s1hat);
        t1.reduce();
        t1.invnttTomont();

        // Add error vector s2
        t1.add(&t1, &s2);

        // Extract t1 and t0
        t1.caddQ();
        var t0: PolyVecK = .{};
        PolyVecK.power2Round(&t1, &t0, &t1);

        // Build public key
        var pk = packing.PublicKey{};
        @memcpy(&pk.rho, rho);
        pk.t1 = t1;

        // Compute H(rho, t1) = tr
        var tr: [TRBYTES]u8 = undefined;
        const pk_bytes = pk.toBytes();
        var tr_hasher = Shake256.init(.{});
        tr_hasher.update(&pk_bytes);
        tr_hasher.squeeze(&tr);

        // Build secret key
        var sk = packing.SecretKey{};
        @memcpy(&sk.rho, rho);
        @memcpy(&sk.key, key);
        @memcpy(&sk.tr, &tr);
        sk.s1 = s1;
        sk.s2 = s2;
        sk.t0 = t0;

        return .{
            .public_key = pk,
            .secret_key = sk,
        };
    }

    pub fn fromSecretKey(sk: packing.SecretKey) Self {
        const rho = &sk.rho;

        // Expand matrix
        var mat: Mat = .{};
        mat.expand(rho);

        // Matrix-vector multiplication
        var s1hat = sk.s1;
        s1hat.ntt();

        var t1: PolyVecK = .{};
        Mat.pointwiseMontgomery(&t1, &mat, &s1hat);
        t1.reduce();
        t1.invnttTomont();

        // Add error vector s2
        t1.add(&t1, &sk.s2);

        // Extract t1 and t0
        t1.caddQ();
        var t0: PolyVecK = .{};
        PolyVecK.power2Round(&t1, &t0, &t1);

        // Build public key
        var pk = packing.PublicKey{};
        @memcpy(&pk.rho, rho);
        pk.t1 = t1;

        return .{
            .public_key = pk,
            .secret_key = sk,
        };
    }

    /// Generate from packed bytes.
    pub fn fromBytes(
        pk_bytes: *const [CRYPTO_PUBLICKEYBYTES]u8,
        sk_bytes: *const [CRYPTO_SECRETKEYBYTES]u8,
    ) Self {
        return .{
            .public_key = packing.PublicKey.fromBytes(pk_bytes),
            .secret_key = packing.SecretKey.fromBytes(sk_bytes),
        };
    }
};

/// Sign a message with optional context.
pub fn sign(
    sk: *const packing.SecretKey,
    msg: []const u8,
    ctx: []const u8,
    random: ?*const [RNDBYTES]u8,
) SigningError![CRYPTO_BYTES]u8 {
    // Validate inputs
    if (msg.len > 1 << 30) return error.MessageTooLong;
    if (ctx.len > 255) return error.ContextTooLong;

    // Prepare prefix = (0, ctxlen, ctx)
    var pre: [257]u8 = undefined;
    pre[0] = 0;
    pre[1] = @intCast(ctx.len);
    @memcpy(pre[2..][0..ctx.len], ctx);

    // Get randomness
    var rnd: [RNDBYTES]u8 = undefined;
    if (random) |r| {
        @memcpy(&rnd, r);
    } else if (signing_mode == .randomized) {
        std.crypto.random.bytes(&rnd);
    } else {
        @memset(&rnd, 0);
    }
    if (simd.has_simd) {
        return signInternalOptimized(sk, msg, pre[0 .. 2 + ctx.len], &rnd);
    }
    return signInternal(sk, msg, pre[0 .. 2 + ctx.len], &rnd);
}

/// Internal signing function.
fn signInternal(
    sk: *const packing.SecretKey,
    msg: []const u8,
    pre: []const u8,
    rnd: *const [RNDBYTES]u8,
) [CRYPTO_BYTES]u8 {
    // Expand matrix
    var mat: Mat = .{};
    mat.expand(&sk.rho);

    // Copy and transform secret vectors
    var s1 = sk.s1;
    var s2 = sk.s2;
    var t0 = sk.t0;
    defer {
        // Zero sensitive secret vectors after use
        secureClearBytes(std.mem.asBytes(&s1));
        secureClearBytes(std.mem.asBytes(&s2));
        secureClearBytes(std.mem.asBytes(&t0));
    }
    s1.ntt();
    s2.ntt();
    t0.ntt();

    // Compute mu = CRH(tr, pre, msg)
    var mu: [CRHBYTES]u8 = undefined;
    defer secureClearBytes(&mu); // Zero hash of secret tr
    {
        var hasher = Shake256.init(.{});
        hasher.update(&sk.tr);
        hasher.update(pre);
        hasher.update(msg);
        hasher.squeeze(&mu);
    }

    // Compute rhoprime = CRH(key, rnd, mu)
    var rhoprime: [CRHBYTES]u8 = undefined;
    defer secureClearBytes(&rhoprime); // Zero derived secret
    {
        var hasher = Shake256.init(.{});
        hasher.update(&sk.key);
        hasher.update(rnd);
        hasher.update(&mu);
        hasher.squeeze(&rhoprime);
    }

    var nonce: u16 = 0;
    var sig: [CRYPTO_BYTES]u8 = undefined;

    // Rejection sampling loop
    while (true) {
        // Sample intermediate vector y
        var y: PolyVecL = .{};
        y.uniformGamma1(&rhoprime, nonce);
        nonce += 1;

        // Matrix-vector multiplication
        var z = y;
        z.ntt();

        var w1: PolyVecK = .{};
        Mat.pointwiseMontgomery(&w1, &mat, &z);
        w1.reduce();
        w1.invnttTomont();

        // Decompose w and call the random oracle
        w1.caddQ();
        var w0: PolyVecK = .{};
        PolyVecK.decompose(&w1, &w0, &w1);

        // In signInternal function:
        var w1_packed: [@as(usize, K) * POLYW1_PACKEDBYTES]u8 = undefined;
        w1.packW1(&w1_packed);

        // Compute challenge hash
        var c_tilde: [CTILDEBYTES]u8 = undefined;
        {
            var hasher = Shake256.init(.{});
            hasher.update(&mu);
            hasher.update(&w1_packed);
            hasher.squeeze(&c_tilde);
        }

        // Sample challenge polynomial
        var cp: Poly = .{};
        cp.challenge(&c_tilde);
        cp.ntt();

        // Compute z, reject if it reveals secret
        z.pointwisePolyMontgomery(&cp, &s1);
        z.invnttTomont();
        z.add(&z, &y);
        z.reduce();

        if (!z.chknorm(GAMMA1 - BETA)) continue;

        // Check that subtracting cs2 does not change high bits
        var h: PolyVecK = .{};
        h.pointwisePolyMontgomery(&cp, &s2);
        h.invnttTomont();
        w0.sub(&w0, &h);
        w0.reduce();

        if (!w0.chknorm(GAMMA2 - BETA)) continue;

        // Compute hints for w1
        h.pointwisePolyMontgomery(&cp, &t0);
        h.invnttTomont();
        h.reduce();

        if (!h.chknorm(GAMMA2)) continue;

        w0.add(&w0, &h);
        const n = PolyVecK.makeHint(&h, &w0, &w1);

        if (n > OMEGA) continue;

        // Pack signature
        packing.packSig(&sig, &c_tilde, &z, &h);
        break;
    }

    return sig;
}

// --- Data Structures for Fast Convolution ---

const SparsePolyIndices = struct { pos: [64]u16, pos_len: usize, neg: [64]u16, neg_len: usize };

const ExtendedPoly = struct {
    data: [2 * N]i32,

    /// SIMD-optimized initialization of the [-s, s] buffer.
    fn init(s: *const Poly) ExtendedPoly {
        var self: ExtendedPoly = undefined;
        // We process in chunks of 8 to use SIMD stores/negation.
        var i: usize = 0;
        while (i < N) : (i += 8) {
            const v = simd.load(@ptrCast(&s.coeffs[i]));
            // Store -s at index i
            const neg_v = simd.broadcast(0) - v;
            simd.store(@ptrCast(&self.data[i]), neg_v);
            // Store s at index N+i
            simd.store(@ptrCast(&self.data[N + i]), v);
        }
        return self;
    }
};

const ExtendedPolyVecK = struct {
    vec: [K]ExtendedPoly,

    fn init(s: *const PolyVecK) ExtendedPolyVecK {
        var self: ExtendedPolyVecK = undefined;
        var k: usize = 0;
        while (k < K) : (k += 1) {
            self.vec[k] = ExtendedPoly.init(&s.vec[k]);
        }
        return self;
    }
};

const ExtendedPolyVecL = struct {
    vec: [L]ExtendedPoly,

    fn init(s: *const PolyVecL) ExtendedPolyVecL {
        var self: ExtendedPolyVecL = undefined;
        var l: usize = 0;
        while (l < L) : (l += 1) {
            self.vec[l] = ExtendedPoly.init(&s.vec[l]);
        }
        return self;
    }
};

// Precomputed offset table for sparse polynomial multiplication
const SparseOffset = struct {
    offset: u16, // Precomputed: N + chunk_start - position
    sign: i32, // +1 for pos, -1 for neg
};

const SparseOffsets = struct {
    offsets: [512]SparseOffset, // Increased size to handle worst-case scenarios
    offsets_len: usize,
};

/// Parse sparse indices once into a unified offset table
/// This eliminates runtime bounds checking in hot loops
/// Precomputes base offsets: base_offset = N - position for each sparse coefficient
inline fn parseSparseOffsets(c: *const Poly, chunk_size: usize, offsets: *SparseOffsets) void {
    offsets.offsets_len = 0;
    for (c.coeffs, 0..) |coeff, i| {
        if (coeff == 1 or coeff == -1) {
            const base_offset = N - i; // Precompute N - position
            if (base_offset + chunk_size <= 2 * N and offsets.offsets_len < offsets.offsets.len) {
                offsets.offsets[offsets.offsets_len] = .{
                    .offset = @intCast(base_offset),
                    .sign = if (coeff == 1) @as(i32, 1) else @as(i32, -1),
                };
                offsets.offsets_len += 1;
            }
        }
    }
}

// --- Optimized Signing Function ---

fn signInternalOptimized(
    sk: *const packing.SecretKey,
    msg: []const u8,
    pre: []const u8,
    rnd: *const [RNDBYTES]u8,
) [CRYPTO_BYTES]u8 {
    // Expand matrix
    var mat: Mat = .{};
    mat.expand(&sk.rho);

    // 1. PRE-COMPUTATION (SIMD Optimized)
    const s1_expanded = ExtendedPolyVecL.init(&sk.s1);
    const s2_expanded = ExtendedPolyVecK.init(&sk.s2);

    // t0 is needed in NTT form for the final Hint check.
    var t0_ntt = sk.t0;
    t0_ntt.ntt();

    // Hashing setup
    var mu: [CRHBYTES]u8 = undefined;
    defer secureClearBytes(&mu);
    {
        var hasher = Shake256.init(.{});
        hasher.update(&sk.tr);
        hasher.update(pre);
        hasher.update(msg);
        hasher.squeeze(&mu);
    }

    var c_tilde_hasher_state = Shake256.init(.{});
    c_tilde_hasher_state.update(&mu);

    var rhoprime: [CRHBYTES]u8 = undefined;
    defer secureClearBytes(&rhoprime);
    {
        var hasher = Shake256.init(.{});
        hasher.update(&sk.key);
        hasher.update(rnd);
        hasher.update(&mu);
        hasher.squeeze(&rhoprime);
    }

    var nonce: u16 = 0;
    var sig: [CRYPTO_BYTES]u8 = undefined;

    // Output buffers
    var z_final: PolyVecL = undefined;
    var h_final: PolyVecK = undefined;
    var c_tilde_final: [CTILDEBYTES]u8 = undefined;
    var c_indices: SparsePolyIndices = undefined;

    // Rejection sampling loop
    while (true) {
        // Sample intermediate vector y
        var y: PolyVecL = undefined;
        y.uniformGamma1(&rhoprime, nonce);
        nonce += 1;

        // Matrix-vector multiplication (A * y)
        var z = y;
        z.ntt();

        var w1: PolyVecK = undefined;
        Mat.pointwiseMontgomery(&w1, &mat, &z);
        w1.reduce();
        w1.invnttTomont();

        // Decompose w
        w1.caddQ();
        var w0: PolyVecK = undefined;
        PolyVecK.decompose(&w1, &w0, &w1);

        var w1_packed: [@as(usize, K) * POLYW1_PACKEDBYTES]u8 = undefined;
        w1.packW1(&w1_packed);

        // Challenge Generation
        var c_tilde: [CTILDEBYTES]u8 = undefined;
        {
            var hasher = c_tilde_hasher_state;
            hasher.update(&w1_packed);
            hasher.squeeze(&c_tilde);
        }

        var cp: Poly = undefined;
        cp.challenge(&c_tilde);

        // Parse sparse coefficients
        parseSparseIndices(&cp, &c_indices);

        // Precompute offset table for optimized PSPM functions
        var c_offsets: SparseOffsets = undefined;
        parseSparseOffsets(&cp, 32, &c_offsets);

        // ---------------------------------------------------------
        // PSPM-TEE: Optimized Checks with Precomputed Offsets
        // ---------------------------------------------------------

        // 3. CHECK z = y + cs1
        var z_valid: bool = undefined;
        if (config.enable_software_pipelining) {
            z_valid = pspmCheckNormSumPipelined(&c_offsets, &s1_expanded, &y, GAMMA1 - BETA, &z_final);
        } else {
            z_valid = pspmCheckNormSumOptimized(&c_offsets, &s1_expanded, &y, GAMMA1 - BETA, &z_final);
        }
        if (!z_valid) {
            continue;
        }

        // 2. CHECK r0 = LowBits(w - cs2)
        var r0_valid: bool = undefined;
        var r_val_cache: PolyVecK = undefined;
        if (config.enable_software_pipelining) {
            r0_valid = pspmCheckNormDiffPipelined(&c_offsets, &s2_expanded, &w0, GAMMA2 - BETA, &r_val_cache);
        } else {
            r0_valid = pspmCheckNormDiffOptimized(&c_offsets, &s2_expanded, &w0, GAMMA2 - BETA, &r_val_cache);
        }
        if (!r0_valid) {
            continue;
        }

        // ---------------------------------------------------------
        // Final Hint Check
        // ---------------------------------------------------------
        cp.ntt();

        var ct0: PolyVecK = undefined;
        ct0.pointwisePolyMontgomery(&cp, &t0_ntt);
        ct0.invnttTomont();
        ct0.reduce();

        if (!ct0.chknorm(GAMMA2)) continue;

        // Recompute rigorous r using optimized diff with precomputed offsets
        var r_val: PolyVecK = r_val_cache;
        r_val.add(&r_val, &ct0);

        const n = PolyVecK.makeHint(&h_final, &r_val, &w1);
        if (n > OMEGA) continue;

        c_tilde_final = c_tilde;
        break;
    }

    packing.packSig(&sig, &c_tilde_final, &z_final, &h_final);
    return sig;
}

// ------------------------------------------------------------------------
// UNROLLED SIMD PSPM Functions
// ------------------------------------------------------------------------

inline fn parseSparseIndices(c: *const Poly, indices: *SparsePolyIndices) void {
    indices.pos_len = 0;
    indices.neg_len = 0;
    for (c.coeffs, 0..) |coeff, i| {
        const is_pos = @intFromBool(coeff == 1);
        const is_neg = @intFromBool(coeff == -1);

        indices.pos[indices.pos_len] = @intCast(i);
        indices.pos_len += is_pos;

        indices.neg[indices.neg_len] = @intCast(i);
        indices.neg_len += is_neg;
    }
}

/// PSPM Check r0.
/// Unrolled 4x (32 coefficients per iteration).
fn pspmCheckNormDiffUnrolled(c: *const SparsePolyIndices, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, limit: i32, dest: *PolyVecK) bool {
    var k: usize = 0;
    var failed: u8 = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        // Iterate in chunks of 32 coefficients (4 SIMD registers)
        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Bounds check before SIMD operations
            if (i + 32 > N) break;

            // Load w0 chunks with bounds checking
            var acc0: simd.I32x8 = undefined;
            var acc1: simd.I32x8 = undefined;
            var acc2: simd.I32x8 = undefined;
            var acc3: simd.I32x8 = undefined;

            if (i + 8 <= N) {
                acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            } else {
                // Fallback for partial chunk
                var temp: [8]i32 = undefined;
                @memset(&temp, 0);
                for (i + 0..@min(i + 8, N)) |j| {
                    temp[j - (i + 0)] = w0_coeffs[j];
                }
                acc0 = simd.load(@ptrCast(&temp));
            }

            if (i + 16 <= N) {
                acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            } else {
                var temp: [8]i32 = undefined;
                @memset(&temp, 0);
                for (i + 8..@min(i + 16, N)) |j| {
                    temp[j - (i + 8)] = w0_coeffs[j];
                }
                acc1 = simd.load(@ptrCast(&temp));
            }

            if (i + 24 <= N) {
                acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            } else {
                var temp: [8]i32 = undefined;
                @memset(&temp, 0);
                for (i + 16..@min(i + 24, N)) |j| {
                    temp[j - (i + 16)] = w0_coeffs[j];
                }
                acc2 = simd.load(@ptrCast(&temp));
            }

            if (i + 32 <= N) {
                acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));
            } else {
                var temp: [8]i32 = undefined;
                @memset(&temp, 0);
                for (i + 24..@min(i + 32, N)) |j| {
                    temp[j - (i + 24)] = w0_coeffs[j];
                }
                acc3 = simd.load(@ptrCast(&temp));
            }

            // Sparse Accumulation Loop with bounds checking
            for (c.pos[0..c.pos_len]) |p| {
                const offset = N + i - p;
                // Bounds check for s_data access
                if (offset + 32 > 2 * N) continue;

                if (offset + 8 <= 2 * N) {
                    acc0 = acc0 - simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                }
                if (offset + 16 <= 2 * N) {
                    acc1 = acc1 - simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                }
                if (offset + 24 <= 2 * N) {
                    acc2 = acc2 - simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                }
                if (offset + 32 <= 2 * N) {
                    acc3 = acc3 - simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
                }
            }

            for (c.neg[0..c.neg_len]) |p| {
                const offset = N + i - p;
                // Bounds check for s_data access
                if (offset + 32 > 2 * N) continue;

                if (offset + 8 <= 2 * N) {
                    acc0 = acc0 + simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                }
                if (offset + 16 <= 2 * N) {
                    acc1 = acc1 + simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                }
                if (offset + 24 <= 2 * N) {
                    acc2 = acc2 + simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                }
                if (offset + 32 <= 2 * N) {
                    acc3 = acc3 + simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
                }
            }

            // Check results
            failed |= @intFromBool(simd.anyGe(simd.abs(acc0), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc1), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc2), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc3), limit));

            // Store results
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
    return failed == 0;
}

/// PSPM Check z.
/// Unrolled 4x.
fn pspmCheckNormSumUnrolled(c: *const SparsePolyIndices, s_expanded: *const ExtendedPolyVecL, y: *const PolyVecL, limit: i32, z_out: *PolyVecL) bool {
    var l: usize = 0;
    var failed: u8 = 0;
    while (l < L) : (l += 1) {
        const s_data = &s_expanded.vec[l].data;
        const y_coeffs = &y.vec[l].coeffs;
        var z_coeffs = &z_out.vec[l].coeffs;

        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load y chunks - no bounds checking needed since N is divisible by 32
            var acc0 = simd.load(@ptrCast(y_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(y_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(y_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(y_coeffs[i + 24 ..].ptr));

            // Sparse Accumulation Loop - optimized for performance
            for (c.pos[0..c.pos_len]) |p| {
                const offset = N + i - p;
                // Since N is divisible by 32 and p <= TAU, offset is safe
                acc0 = acc0 + simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                acc1 = acc1 + simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                acc2 = acc2 + simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                acc3 = acc3 + simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
            }

            for (c.neg[0..c.neg_len]) |p| {
                const offset = N + i - p;
                // Since N is divisible by 32 and p <= TAU, offset is safe
                acc0 = acc0 - simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                acc1 = acc1 - simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                acc2 = acc2 - simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                acc3 = acc3 - simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
            }

            failed |= @intFromBool(simd.anyGe(simd.abs(acc0), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc1), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc2), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc3), limit));

            // Store results - no bounds checking needed since N is divisible by 32
            simd.store(@ptrCast(z_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(z_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(z_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(z_coeffs[i + 24 ..].ptr), acc3);
        }
    }
    return failed == 0;
}

fn pspmComputeDiffUnrolled(c: *const SparsePolyIndices, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, dest: *PolyVecK) void {
    var k: usize = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load w0 chunks - no bounds checking needed since N is divisible by 32
            var acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));

            // Sparse Accumulation Loop - optimized for performance
            for (c.pos[0..c.pos_len]) |p| {
                const offset = N + i - p;
                // Since N is divisible by 32 and p <= TAU, offset is safe
                acc0 = acc0 - simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                acc1 = acc1 - simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                acc2 = acc2 - simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                acc3 = acc3 - simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
            }

            for (c.neg[0..c.neg_len]) |p| {
                const offset = N + i - p;
                // Since N is divisible by 32 and p <= TAU, offset is safe
                acc0 = acc0 + simd.loadU(@ptrCast(s_data[offset + 0 ..].ptr));
                acc1 = acc1 + simd.loadU(@ptrCast(s_data[offset + 8 ..].ptr));
                acc2 = acc2 + simd.loadU(@ptrCast(s_data[offset + 16 ..].ptr));
                acc3 = acc3 + simd.loadU(@ptrCast(s_data[offset + 24 ..].ptr));
            }

            // Store results - no bounds checking needed since N is divisible by 32
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
}

// ------------------------------------------------------------------------
// OPTIMIZED PSPM FUNCTIONS WITH PRECOMPUTED OFFSET TABLE
// ------------------------------------------------------------------------

/// PSPM Check r0 using precomputed offset table
/// Eliminates runtime bounds checking in hot loops
fn pspmCheckNormDiffOptimized(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, limit: i32, dest: *PolyVecK) bool {
    var k: usize = 0;
    var failed: u8 = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        // Iterate in chunks of 32 coefficients (4 SIMD registers)
        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load w0 chunks - no bounds checking needed since N is divisible by 32
            var acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));

            // Optimized sparse accumulation using precomputed offsets
            for (offsets.offsets[0..offsets.offsets_len]) |offset| {
                const offset_idx = offset.offset + i; // base_offset + i = (N - position) + i = N + i - position
                const sign_vec = simd.broadcast(offset.sign);

                // Bounds check to ensure we don't access beyond s_data array
                if (offset_idx + 32 <= 2 * N) {
                    acc0 = acc0 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    acc1 = acc1 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                    acc2 = acc2 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                    acc3 = acc3 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));
                }
            }

            failed |= @intFromBool(simd.anyGe(simd.abs(acc0), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc1), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc2), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc3), limit));

            // Store results - no bounds checking needed since N is divisible by 32
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
    return failed == 0;
}

/// PSPM Check z using precomputed offset table
fn pspmCheckNormSumOptimized(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecL, y: *const PolyVecL, limit: i32, z_out: *PolyVecL) bool {
    const limit_vec = simd.broadcast(@as(i32, limit));
    const sign_pos = simd.broadcast(@as(i32, 1));
    const sign_neg = simd.broadcast(@as(i32, -1));

    const offsets_slice = offsets.offsets[0..offsets.offsets_len];

    var failed: u8 = 0;
    var l: usize = 0;

    while (l < L) : (l += 1) {
        const s_data = &s_expanded.vec[l].data;
        const y_coeffs = &y.vec[l].coeffs;
        var z_coeffs = &z_out.vec[l].coeffs;

        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load y chunks - no bounds checking needed since N is divisible by 32
            var acc0 = simd.load(@ptrCast(y_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(y_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(y_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(y_coeffs[i + 24 ..].ptr));
            const max_valid_offset = 2 * N - 32 - i;
            // Optimized sparse accumulation using precomputed offsets
            for (offsets_slice) |offset| {
                if (offset.offset > max_valid_offset) {
                    break;
                }
                const offset_idx = offset.offset + i; // base_offset + i = (N - position) + i = N + i - position
                const sign_vec = if (offset.sign == 1) sign_pos else sign_neg;
                acc0 += sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                acc1 += sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                acc2 += sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                acc3 += sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));
            }

            failed |= @intFromBool(simd.anyGeVec(simd.abs(acc0), limit_vec));
            failed |= @intFromBool(simd.anyGeVec(simd.abs(acc1), limit_vec));
            failed |= @intFromBool(simd.anyGeVec(simd.abs(acc2), limit_vec));
            failed |= @intFromBool(simd.anyGeVec(simd.abs(acc3), limit_vec));

            // Store results - no bounds checking needed since N is divisible by 32
            simd.store(@ptrCast(z_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(z_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(z_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(z_coeffs[i + 24 ..].ptr), acc3);
        }
    }
    return failed == 0;
}

/// PSPM Compute diff using precomputed offset table
fn pspmComputeDiffOptimized(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, dest: *PolyVecK) void {
    var k: usize = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        var i: usize = 0;
        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load w0 chunks - no bounds checking needed since N is divisible by 32
            var acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));

            // Optimized sparse accumulation using precomputed offsets
            for (offsets.offsets[0..offsets.offsets_len]) |offset| {
                const offset_idx = offset.offset + i; // base_offset + i = (N - position) + i = N + i - position
                const sign_vec = simd.broadcast(offset.sign);

                // Bounds check to ensure we don't access beyond s_data array
                if (offset_idx + 32 <= 2 * N) {
                    acc0 = acc0 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    acc1 = acc1 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                    acc2 = acc2 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                    acc3 = acc3 - sign_vec * simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));
                }
            }

            // Store results - no bounds checking needed since N is divisible by 32
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
}

/// Verify a signature with optional context.
pub fn verify(
    pk: *const packing.PublicKey,
    msg: []const u8,
    ctx: []const u8,
    sig: *const [CRYPTO_BYTES]u8,
) SigningError!void {
    // Validate inputs
    if (msg.len > 1 << 30) return error.MessageTooLong;
    if (ctx.len > 255) return error.ContextTooLong;

    // Prepare prefix = (0, ctxlen, ctx)
    var pre: [257]u8 = undefined;
    pre[0] = 0;
    pre[1] = @intCast(ctx.len);
    @memcpy(pre[2..][0..ctx.len], ctx);

    return verifyInternal(pk, msg, pre[0 .. 2 + ctx.len], sig);
}

/// Internal verification function (constant-time).
fn verifyInternal(
    pk: *const packing.PublicKey,
    msg: []const u8,
    pre: []const u8,
    sig: *const [CRYPTO_BYTES]u8,
) SigningError!void {
    // Unpack signature - accumulate failure status instead of early return
    var c: [CTILDEBYTES]u8 = undefined;
    var z: PolyVecL = .{};
    var h: PolyVecK = .{};

    const unpack_ok = packing.unpackSig(&c, &z, &h, sig);
    var valid = @intFromBool(unpack_ok);

    // Check z norm - accumulate failure instead of early return
    const z_valid = z.chknorm(GAMMA1 - BETA);
    valid &= @intFromBool(z_valid);

    // Compute CRH(H(rho, t1), pre, msg) - always compute regardless of previous failures
    var mu: [CRHBYTES]u8 = undefined;
    {
        const pk_bytes = pk.toBytes();
        var tr: [TRBYTES]u8 = undefined;
        var hasher = Shake256.init(.{});
        hasher.update(&pk_bytes);
        hasher.squeeze(&tr);

        hasher = Shake256.init(.{});
        hasher.update(&tr);
        hasher.update(pre);
        hasher.update(msg);
        hasher.squeeze(&mu);
    }

    // Sample challenge polynomial - always sample
    var cp: Poly = .{};
    cp.challenge(&c);

    // Expand matrix - always expand
    var mat: Mat = .{};
    mat.expand(&pk.rho);

    // Compute Az - c*2^d*t1 - always compute full operation
    z.ntt();

    var w1: PolyVecK = .{};
    Mat.pointwiseMontgomery(&w1, &mat, &z);

    cp.ntt();

    var t1 = pk.t1;
    t1.shiftl();
    t1.ntt();
    t1.pointwisePolyMontgomery(&cp, &t1);

    w1.sub(&w1, &t1);
    w1.reduce();
    w1.invnttTomont();

    // Reconstruct w1 - always reconstruct
    w1.caddQ();
    PolyVecK.useHint(&w1, &w1, &h);

    // Pack w1 - always pack
    var w1_packed: [@as(usize, K) * POLYW1_PACKEDBYTES]u8 = undefined;
    w1.packW1(&w1_packed);

    // Call random oracle and verify challenge - always compute
    var c2: [CTILDEBYTES]u8 = undefined;
    {
        var hasher = Shake256.init(.{});
        hasher.update(&mu);
        hasher.update(&w1_packed);
        hasher.squeeze(&c2);
    }

    // Final challenge comparison - accumulate result
    const challenge_valid = std.mem.eql(u8, &c, &c2);
    valid &= @intFromBool(challenge_valid);

    // Return based on accumulated validity - single point of failure
    if (valid == 0) {
        return error.VerificationFailed;
    }
}

/// Sign a message and prepend signature to output.
pub fn signMessage(
    sk: *const packing.SecretKey,
    msg: []const u8,
    ctx: []const u8,
    out: []u8,
) SigningError!usize {
    if (ctx.len > 255) return error.ContextTooLong;
    if (out.len < CRYPTO_BYTES + msg.len) return error.InvalidSignatureLength;

    // Copy message to end of output (backwards for overlapping buffers)
    var i: usize = msg.len;
    while (i > 0) {
        i -= 1;
        out[CRYPTO_BYTES + i] = msg[i];
    }

    // Sign
    const sig = try sign(sk, out[CRYPTO_BYTES..][0..msg.len], ctx, null);
    @memcpy(out[0..CRYPTO_BYTES], &sig);

    return CRYPTO_BYTES + msg.len;
}

/// Verify and extract message from signed message.
pub fn openMessage(
    pk: *const packing.PublicKey,
    sm: []const u8,
    ctx: []const u8,
    out: []u8,
) SigningError!usize {
    if (sm.len < CRYPTO_BYTES) return error.InvalidSignatureLength;

    const mlen = sm.len - CRYPTO_BYTES;

    try verify(pk, sm[CRYPTO_BYTES..], ctx, sm[0..CRYPTO_BYTES]);

    // Copy message to output
    @memcpy(out[0..mlen], sm[CRYPTO_BYTES..][0..mlen]);

    return mlen;
}

// ==================== C-Compatible API ====================

/// Generate key pair (C-compatible).
pub fn cryptoSignKeypair(
    pk: *[CRYPTO_PUBLICKEYBYTES]u8,
    sk: *[CRYPTO_SECRETKEYBYTES]u8,
) c_int {
    const kp = KeyPair.generate(null);
    kp.public_key.pack(pk);
    kp.secret_key.pack(sk);
    return 0;
}

/// Sign message (C-compatible).
pub fn cryptoSignSignature(
    sig: *[CRYPTO_BYTES]u8,
    siglen: *usize,
    m: [*]const u8,
    mlen: usize,
    ctx: [*]const u8,
    ctxlen: usize,
    sk: *const [CRYPTO_SECRETKEYBYTES]u8,
) c_int {
    if (ctxlen > 255) return -1;

    const secret_key = packing.SecretKey.fromBytes(sk);
    const result = sign(&secret_key, m[0..mlen], ctx[0..ctxlen], null) catch return -1;

    @memcpy(sig, &result);
    siglen.* = CRYPTO_BYTES;
    return 0;
}

/// Verify signature (C-compatible).
pub fn cryptoSignVerify(
    sig: *const [CRYPTO_BYTES]u8,
    siglen: usize,
    m: [*]const u8,
    mlen: usize,
    ctx: [*]const u8,
    ctxlen: usize,
    pk: *const [CRYPTO_PUBLICKEYBYTES]u8,
) c_int {
    if (siglen != CRYPTO_BYTES) return -1;
    if (ctxlen > 255) return -1;

    const public_key = packing.PublicKey.fromBytes(pk);
    verify(&public_key, m[0..mlen], ctx[0..ctxlen], sig) catch return -1;
    return 0;
}

/// Sign message with prepended signature (C-compatible).
pub fn cryptoSign(
    sm: [*]u8,
    smlen: *usize,
    m: [*]const u8,
    mlen: usize,
    ctx: [*]const u8,
    ctxlen: usize,
    sk: *const [CRYPTO_SECRETKEYBYTES]u8,
) c_int {
    if (ctxlen > 255) return -1;

    const secret_key = packing.SecretKey.fromBytes(sk);

    // Copy message backwards
    var i: usize = mlen;
    while (i > 0) {
        i -= 1;
        sm[CRYPTO_BYTES + i] = m[i];
    }

    const sig = sign(&secret_key, sm[CRYPTO_BYTES..][0..mlen], ctx[0..ctxlen], null) catch return -1;
    @memcpy(sm[0..CRYPTO_BYTES], &sig);
    smlen.* = CRYPTO_BYTES + mlen;
    return 0;
}

/// Open signed message (C-compatible).
pub fn cryptoSignOpen(
    m: [*]u8,
    mlen: *usize,
    sm: [*]const u8,
    smlen: usize,
    ctx: [*]const u8,
    ctxlen: usize,
    pk: *const [CRYPTO_PUBLICKEYBYTES]u8,
) c_int {
    if (smlen < CRYPTO_BYTES) {
        mlen.* = 0;
        return -1;
    }
    if (ctxlen > 255) {
        mlen.* = 0;
        return -1;
    }

    const public_key = packing.PublicKey.fromBytes(pk);
    const msg_len = smlen - CRYPTO_BYTES;

    verify(&public_key, sm[CRYPTO_BYTES..][0..msg_len], ctx[0..ctxlen], sm[0..CRYPTO_BYTES]) catch {
        mlen.* = 0;
        secureClearBytes(m[0..smlen]);
        return -1;
    };

    @memcpy(m[0..msg_len], sm[CRYPTO_BYTES..][0..msg_len]);
    mlen.* = msg_len;
    return 0;
}

// ------------------------------------------------------------------------
// SOFTWARE PIPELINED PSPM FUNCTIONS - HIDE LOAD LATENCY
// ------------------------------------------------------------------------

/// Software pipelined PSPM Check r0 - hides load latency through prefetching
/// Uses double-buffering technique to overlap computation with memory access
fn pspmCheckNormDiffPipelined(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, limit: i32, dest: *PolyVecK) bool {
    var k: usize = 0;
    var failed: u8 = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        // Iterate in chunks of 32 coefficients (4 SIMD registers)
        var i: usize = 0;

        // Prologue: Prefetch first iteration data
        if (N >= 32) {
            // Prefetch w0 data for first iteration
            _ = simd.load(@ptrCast(w0_coeffs[0..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[8..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[16..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[24..].ptr));

            // Prefetch s_data for common offsets (first few offsets are typically hot)
            const prefetch_offsets = @min(offsets.offsets_len, 4);
            for (offsets.offsets[0..prefetch_offsets]) |offset| {
                const offset_idx = offset.offset;
                if (offset_idx + 32 <= 2 * N) {
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                }
            }
        }

        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load w0 chunks with prefetch for next iteration
            var acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));

            // Prefetch next iteration w0 data (if available)
            if (i + 64 <= N) {
                _ = simd.load(@ptrCast(w0_coeffs[i + 32 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 40 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 48 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 56 ..].ptr));
            }

            // Optimized sparse accumulation with software pipelining
            for (offsets.offsets[0..offsets.offsets_len]) |offset| {
                const offset_idx = offset.offset + i; // base_offset + i = (N - position) + i = N + i - position
                const sign_vec = simd.broadcast(offset.sign);

                // Bounds check to ensure we don't access beyond s_data array
                if (offset_idx + 32 <= 2 * N) {
                    // Overlap: Start loading s_data while previous computation is in flight
                    const s_chunk0 = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    const s_chunk1 = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                    const s_chunk2 = simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                    const s_chunk3 = simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));

                    // Computation happens while next data is being fetched
                    acc0 = acc0 - sign_vec * s_chunk0;
                    acc1 = acc1 - sign_vec * s_chunk1;
                    acc2 = acc2 - sign_vec * s_chunk2;
                    acc3 = acc3 - sign_vec * s_chunk3;
                }
            }

            // Check norms - this creates a small bubble but is necessary
            failed |= @intFromBool(simd.anyGe(simd.abs(acc0), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc1), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc2), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc3), limit));

            // Store results
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
    return failed == 0;
}

/// Software pipelined PSPM Check z - hides load latency through prefetching
fn pspmCheckNormSumPipelined(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecL, y: *const PolyVecL, limit: i32, z_out: *PolyVecL) bool {
    var l: usize = 0;
    var failed: u8 = 0;
    while (l < L) : (l += 1) {
        const s_data = &s_expanded.vec[l].data;
        const y_coeffs = &y.vec[l].coeffs;
        var z_coeffs = &z_out.vec[l].coeffs;

        var i: usize = 0;

        // Prologue: Prefetch first iteration data
        if (N >= 32) {
            _ = simd.load(@ptrCast(y_coeffs[0..].ptr));
            _ = simd.load(@ptrCast(y_coeffs[8..].ptr));
            _ = simd.load(@ptrCast(y_coeffs[16..].ptr));
            _ = simd.load(@ptrCast(y_coeffs[24..].ptr));

            const prefetch_offsets = @min(offsets.offsets_len, 4);
            for (offsets.offsets[0..prefetch_offsets]) |offset| {
                const offset_idx = offset.offset;
                if (offset_idx + 32 <= 2 * N) {
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                }
            }
        }

        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load y chunks with prefetch for next iteration
            var acc0 = simd.load(@ptrCast(y_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(y_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(y_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(y_coeffs[i + 24 ..].ptr));

            // Prefetch next iteration y data
            if (i + 64 <= N) {
                _ = simd.load(@ptrCast(y_coeffs[i + 32 ..].ptr));
                _ = simd.load(@ptrCast(y_coeffs[i + 40 ..].ptr));
                _ = simd.load(@ptrCast(y_coeffs[i + 48 ..].ptr));
                _ = simd.load(@ptrCast(y_coeffs[i + 56 ..].ptr));
            }

            // Optimized sparse accumulation with software pipelining
            for (offsets.offsets[0..offsets.offsets_len]) |offset| {
                const offset_idx = offset.offset + i;
                const sign_vec = simd.broadcast(offset.sign);

                if (offset_idx + 32 <= 2 * N) {
                    // Overlap: Load s_data while previous computation is in flight
                    const s_chunk0 = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    const s_chunk1 = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                    const s_chunk2 = simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                    const s_chunk3 = simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));

                    // Computation happens while next data is being fetched
                    acc0 = acc0 + sign_vec * s_chunk0;
                    acc1 = acc1 + sign_vec * s_chunk1;
                    acc2 = acc2 + sign_vec * s_chunk2;
                    acc3 = acc3 + sign_vec * s_chunk3;
                }
            }

            // Store results with pipelining - store current while prefetching next
            simd.store(@ptrCast(z_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(z_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(z_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(z_coeffs[i + 24 ..].ptr), acc3);

            // Check norms
            failed |= @intFromBool(simd.anyGe(simd.abs(acc0), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc1), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc2), limit));
            failed |= @intFromBool(simd.anyGe(simd.abs(acc3), limit));
        }
    }
    return failed == 0;
}

/// Software pipelined PSPM compute difference - hides load latency
fn pspmComputeDiffPipelined(offsets: *const SparseOffsets, s_expanded: *const ExtendedPolyVecK, w0: *const PolyVecK, dest: *PolyVecK) void {
    var k: usize = 0;
    while (k < K) : (k += 1) {
        const s_data = &s_expanded.vec[k].data;
        const w0_coeffs = &w0.vec[k].coeffs;
        var d_coeffs = &dest.vec[k].coeffs;

        var i: usize = 0;

        // Prologue: Prefetch first iteration data
        if (N >= 32) {
            _ = simd.load(@ptrCast(w0_coeffs[0..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[8..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[16..].ptr));
            _ = simd.load(@ptrCast(w0_coeffs[24..].ptr));

            const prefetch_offsets = @min(offsets.offsets_len, 4);
            for (offsets.offsets[0..prefetch_offsets]) |offset| {
                const offset_idx = offset.offset;
                if (offset_idx + 32 <= 2 * N) {
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    _ = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                }
            }
        }

        while (i < N) : (i += 32) {
            // Comptime assertion: N is divisible by 32
            comptime assert(N % 32 == 0);

            // Load w0 chunks with prefetch for next iteration
            var acc0 = simd.load(@ptrCast(w0_coeffs[i + 0 ..].ptr));
            var acc1 = simd.load(@ptrCast(w0_coeffs[i + 8 ..].ptr));
            var acc2 = simd.load(@ptrCast(w0_coeffs[i + 16 ..].ptr));
            var acc3 = simd.load(@ptrCast(w0_coeffs[i + 24 ..].ptr));

            // Prefetch next iteration w0 data
            if (i + 64 <= N) {
                _ = simd.load(@ptrCast(w0_coeffs[i + 32 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 40 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 48 ..].ptr));
                _ = simd.load(@ptrCast(w0_coeffs[i + 56 ..].ptr));
            }

            // Optimized sparse accumulation with software pipelining
            for (offsets.offsets[0..offsets.offsets_len]) |offset| {
                const offset_idx = offset.offset + i;
                const sign_vec = simd.broadcast(offset.sign);

                if (offset_idx + 32 <= 2 * N) {
                    // Overlap: Load s_data while previous computation is in flight
                    const s_chunk0 = simd.loadU(@ptrCast(s_data[offset_idx + 0 ..].ptr));
                    const s_chunk1 = simd.loadU(@ptrCast(s_data[offset_idx + 8 ..].ptr));
                    const s_chunk2 = simd.loadU(@ptrCast(s_data[offset_idx + 16 ..].ptr));
                    const s_chunk3 = simd.loadU(@ptrCast(s_data[offset_idx + 24 ..].ptr));

                    // Computation happens while next data is being fetched
                    acc0 = acc0 - sign_vec * s_chunk0;
                    acc1 = acc1 - sign_vec * s_chunk1;
                    acc2 = acc2 - sign_vec * s_chunk2;
                    acc3 = acc3 - sign_vec * s_chunk3;
                }
            }

            // Store results - overlap store with next iteration prefetch
            simd.store(@ptrCast(d_coeffs[i + 0 ..].ptr), acc0);
            simd.store(@ptrCast(d_coeffs[i + 8 ..].ptr), acc1);
            simd.store(@ptrCast(d_coeffs[i + 16 ..].ptr), acc2);
            simd.store(@ptrCast(d_coeffs[i + 24 ..].ptr), acc3);
        }
    }
}

// ==================== Tests ====================

test "keypair generation" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    // Verify we can pack/unpack
    const pk_bytes = kp.public_key.toBytes();
    const sk_bytes = kp.secret_key.toBytes();

    const kp2 = KeyPair.fromBytes(&pk_bytes, &sk_bytes);
    try std.testing.expectEqualSlices(u8, &kp.public_key.rho, &kp2.public_key.rho);
    try std.testing.expectEqualSlices(u8, &kp.secret_key.key, &kp2.secret_key.key);
}

test "sign and verify" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const ctx = "context";

    const sig = try sign(&kp.secret_key, msg, ctx, null);
    try verify(&kp.public_key, msg, ctx, &sig);
}

test "sign and verify empty context" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";

    const sig = try sign(&kp.secret_key, msg, "", null);
    try verify(&kp.public_key, msg, "", &sig);
}

test "verify fails with wrong message" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const ctx = "";

    const sig = try sign(&kp.secret_key, msg, ctx, null);

    const result = verify(&kp.public_key, "wrong message", ctx, &sig);
    try std.testing.expectError(error.VerificationFailed, result);
}

test "verify fails with wrong context" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";

    const sig = try sign(&kp.secret_key, msg, "ctx1", null);

    const result = verify(&kp.public_key, msg, "ctx2", &sig);
    try std.testing.expectError(error.VerificationFailed, result);
}

test "context too long" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const long_ctx = [_]u8{0} ** 256;

    const result = sign(&kp.secret_key, msg, &long_ctx, null);
    try std.testing.expectError(error.ContextTooLong, result);
}

test "sign message and open" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const ctx = "";

    var sm: [CRYPTO_BYTES + msg.len]u8 = undefined;
    const smlen = try signMessage(&kp.secret_key, msg, ctx, &sm);

    var out: [msg.len]u8 = undefined;
    const mlen = try openMessage(&kp.public_key, sm[0..smlen], ctx, &out);

    try std.testing.expectEqualSlices(u8, msg, out[0..mlen]);
}
