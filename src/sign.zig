//! Dilithium signature scheme.
//! ML-DSA (FIPS 204) implementation.

const std = @import("std");
const params = @import("params.zig");
const poly_mod = @import("poly.zig");
const polyvec_mod = @import("polyvec.zig");
const packing = @import("packing.zig");
const config = @import("config.zig");

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
        defer @memset(&seedbuf, 0); // Zero sensitive seed material

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
) SigningError![CRYPTO_BYTES]u8 {
    if (ctx.len > 255) return error.ContextTooLong;

    // Prepare prefix = (0, ctxlen, ctx)
    var pre: [257]u8 = undefined;
    pre[0] = 0;
    pre[1] = @intCast(ctx.len);
    @memcpy(pre[2..][0..ctx.len], ctx);

    // Get randomness
    var rnd: [RNDBYTES]u8 = undefined;
    if (signing_mode == .randomized) {
        std.crypto.random.bytes(&rnd);
    } else {
        @memset(&rnd, 0);
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
        @memset(std.mem.asBytes(&s1), 0);
        @memset(std.mem.asBytes(&s2), 0);
        @memset(std.mem.asBytes(&t0), 0);
    }
    s1.ntt();
    s2.ntt();
    t0.ntt();

    // Compute mu = CRH(tr, pre, msg)
    var mu: [CRHBYTES]u8 = undefined;
    defer @memset(&mu, 0); // Zero hash of secret tr
    {
        var hasher = Shake256.init(.{});
        hasher.update(&sk.tr);
        hasher.update(pre);
        hasher.update(msg);
        hasher.squeeze(&mu);
    }

    // Compute rhoprime = CRH(key, rnd, mu)
    var rhoprime: [CRHBYTES]u8 = undefined;
    defer @memset(&rhoprime, 0); // Zero derived secret
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

/// Verify a signature with optional context.
pub fn verify(
    pk: *const packing.PublicKey,
    msg: []const u8,
    ctx: []const u8,
    sig: *const [CRYPTO_BYTES]u8,
) SigningError!void {
    if (ctx.len > 255) return error.ContextTooLong;

    // Prepare prefix = (0, ctxlen, ctx)
    var pre: [257]u8 = undefined;
    pre[0] = 0;
    pre[1] = @intCast(ctx.len);
    @memcpy(pre[2..][0..ctx.len], ctx);

    return verifyInternal(pk, msg, pre[0 .. 2 + ctx.len], sig);
}

/// Internal verification function.
fn verifyInternal(
    pk: *const packing.PublicKey,
    msg: []const u8,
    pre: []const u8,
    sig: *const [CRYPTO_BYTES]u8,
) SigningError!void {
    // Unpack signature
    var c: [CTILDEBYTES]u8 = undefined;
    var z: PolyVecL = .{};
    var h: PolyVecK = .{};

    if (!packing.unpackSig(&c, &z, &h, sig)) {
        return error.InvalidSignature;
    }

    if (!z.chknorm(GAMMA1 - BETA)) {
        return error.InvalidSignature;
    }

    // Compute CRH(H(rho, t1), pre, msg)
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

    // Sample challenge polynomial
    var cp: Poly = .{};
    cp.challenge(&c);

    // Expand matrix
    var mat: Mat = .{};
    mat.expand(&pk.rho);

    // Compute Az - c*2^d*t1
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

    // Reconstruct w1
    w1.caddQ();
    PolyVecK.useHint(&w1, &w1, &h);

    // In verifyInternal function:
    var w1_packed: [@as(usize, K) * POLYW1_PACKEDBYTES]u8 = undefined;
    w1.packW1(&w1_packed);

    // Call random oracle and verify challenge
    var c2: [CTILDEBYTES]u8 = undefined;
    {
        var hasher = Shake256.init(.{});
        hasher.update(&mu);
        hasher.update(&w1_packed);
        hasher.squeeze(&c2);
    }

    if (!std.mem.eql(u8, &c, &c2)) {
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
    const sig = try sign(sk, out[CRYPTO_BYTES..][0..msg.len], ctx);
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
    const result = sign(&secret_key, m[0..mlen], ctx[0..ctxlen]) catch return -1;

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

    const sig = sign(&secret_key, sm[CRYPTO_BYTES..][0..mlen], ctx[0..ctxlen]) catch return -1;
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
        @memset(m[0..smlen], 0);
        return -1;
    };

    @memcpy(m[0..msg_len], sm[CRYPTO_BYTES..][0..msg_len]);
    mlen.* = msg_len;
    return 0;
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

    const sig = try sign(&kp.secret_key, msg, ctx);
    try verify(&kp.public_key, msg, ctx, &sig);
}

test "sign and verify empty context" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";

    const sig = try sign(&kp.secret_key, msg, "");
    try verify(&kp.public_key, msg, "", &sig);
}

test "verify fails with wrong message" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const ctx = "";

    const sig = try sign(&kp.secret_key, msg, ctx);

    const result = verify(&kp.public_key, "wrong message", ctx, &sig);
    try std.testing.expectError(error.VerificationFailed, result);
}

test "verify fails with wrong context" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";

    const sig = try sign(&kp.secret_key, msg, "ctx1");

    const result = verify(&kp.public_key, msg, "ctx2", &sig);
    try std.testing.expectError(error.VerificationFailed, result);
}

test "context too long" {
    const seed = [_]u8{0x42} ** SEEDBYTES;
    const kp = KeyPair.generate(&seed);

    const msg = "test message";
    const long_ctx = [_]u8{0} ** 256;

    const result = sign(&kp.secret_key, msg, &long_ctx);
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
