//! Symmetric cryptographic primitives for Dilithium.
//! Uses SHAKE128/SHAKE256 for stream generation.

const std = @import("std");
const params = @import("params.zig");

const SEEDBYTES = params.SEEDBYTES;
const CRHBYTES = params.CRHBYTES;

// SHAKE rates (bytes per block)
pub const SHAKE128_RATE = 168;
pub const SHAKE256_RATE = 136;

pub const STREAM128_BLOCKBYTES = SHAKE128_RATE;
pub const STREAM256_BLOCKBYTES = SHAKE256_RATE;

/// Stream state for SHAKE128-based expansion.
pub const Stream128State = struct {
    state: std.crypto.hash.sha3.Shake128,

    const Self = @This();

    /// Initialize the stream with a seed and nonce.
    /// Absorbs: seed || nonce (little-endian)
    pub fn init(seed: *const [SEEDBYTES]u8, nonce: u16) Self {
        var state = std.crypto.hash.sha3.Shake128.init(.{});
        state.update(seed);
        state.update(&std.mem.toBytes(std.mem.nativeToLittle(u16, nonce)));
        return .{ .state = state };
    }

    /// Squeeze a single block of output.
    pub fn squeezeBlock(self: *Self) [SHAKE128_RATE]u8 {
        var block: [SHAKE128_RATE]u8 = undefined;
        self.state.squeeze(&block);
        return block;
    }

    /// Squeeze multiple blocks into output buffer.
    pub fn squeezeBlocks(self: *Self, out: []u8, nblocks: usize) void {
        std.debug.assert(out.len >= nblocks * SHAKE128_RATE);
        const len = nblocks * SHAKE128_RATE;
        self.state.squeeze(out[0..len]);
    }

    /// Squeeze arbitrary number of bytes.
    pub fn squeeze(self: *Self, out: []u8) void {
        self.state.squeeze(out);
    }
};

/// Stream state for SHAKE256-based expansion.
pub const Stream256State = struct {
    state: std.crypto.hash.sha3.Shake256,

    const Self = @This();

    /// Initialize the stream with a seed and nonce.
    /// Absorbs: seed || nonce (little-endian)
    pub fn init(seed: *const [CRHBYTES]u8, nonce: u16) Self {
        var state = std.crypto.hash.sha3.Shake256.init(.{});
        state.update(seed);
        state.update(&std.mem.toBytes(std.mem.nativeToLittle(u16, nonce)));
        return .{ .state = state };
    }

    /// Squeeze a single block of output.
    pub fn squeezeBlock(self: *Self) [SHAKE256_RATE]u8 {
        var block: [SHAKE256_RATE]u8 = undefined;
        self.state.squeeze(&block);
        return block;
    }

    /// Squeeze multiple blocks into output buffer.
    pub fn squeezeBlocks(self: *Self, out: []u8, nblocks: usize) void {
        std.debug.assert(out.len >= nblocks * SHAKE256_RATE);
        const len = nblocks * SHAKE256_RATE;
        self.state.squeeze(out[0..len]);
    }

    /// Squeeze arbitrary number of bytes.
    pub fn squeeze(self: *Self, out: []u8) void {
        self.state.squeeze(out);
    }
};

// Convenience aliases matching C naming
pub const stream128_state = Stream128State;
pub const stream256_state = Stream256State;

/// Initialize SHAKE128 stream (C-compatible naming).
pub fn dilithium_shake128_stream_init(seed: *const [SEEDBYTES]u8, nonce: u16) Stream128State {
    return Stream128State.init(seed, nonce);
}

/// Initialize SHAKE256 stream (C-compatible naming).
pub fn dilithium_shake256_stream_init(seed: *const [CRHBYTES]u8, nonce: u16) Stream256State {
    return Stream256State.init(seed, nonce);
}

const keccak4x = @import("keccak4x.zig");

/// Parallel SHAKE128 state for 4 streams.
pub const Stream128x4State = struct {
    state: keccak4x.State4x,
    offset: usize = 0, // Byte offset in the current block (same for all lanes)

    const Self = @This();

    /// Initialize 4 SHAKE128 states with same seed but distinct nonces.
    pub fn init(seed: *const [SEEDBYTES]u8, nonces: [4]u16) Self {
        var self = Self{
            .state = undefined,
        };

        // Clear state
        @memset(std.mem.asBytes(&self.state.s), 0);

        // Prepare input block: seed || nonce (little-endian)
        // Length = 32 + 2 = 34 bytes.
        // SHAKE128 rate is 168 bytes. Block fits in one absorption.

        // Load seed into first 4 u64 words (32 bytes)
        const seed_u64 = std.mem.bytesAsSlice(u64, seed);
        for (0..4) |i| {
            self.state.s[i] = @splat(seed_u64[i]);
        }

        // Load nonces
        // The seed ended at byte 32.
        // Byte 32 and 33 are the nonce.
        // This corresponds to the beginning of state.s[4] (bytes 32-39).

        // Construct the 5th u64 word for each lane
        var word4: keccak4x.U64x4 = undefined;
        var nonces_arr: [4]u64 = undefined;

        for (0..4) |i| {
            // nonce (16-bit) | 0x1F (8-bit) << 16
            const val: u64 = @as(u64, nonces[i]) | (@as(u64, 0x1F) << 16);
            nonces_arr[i] = val;
        }
        word4 = nonces_arr;

        self.state.s[4] = word4;

        // Add final padding bit (0x80) at the end of the rate (byte 167).
        // Rate is 168 bytes = 21 words of 8 bytes.
        // The last byte is in word 20 (bytes 160-167).
        // Byte 167 is the MSB of word 20.
        // 0x80 << 56

        const pad_final: u64 = @as(u64, 0x80) << 56;
        self.state.s[20] ^= @as(keccak4x.U64x4, @splat(pad_final));

        // Permute
        keccak4x.keccakF1600_4x(&self.state);

        self.offset = 0;
        return self;
    }

    /// Squeeze 4 blocks of output into 4 buffers.
    /// Each buffer must be at least SHAKE128_RATE bytes.
    pub fn squeezeBlocks(self: *Self, out0: []u8, out1: []u8, out2: []u8, out3: []u8) void {
        const rate = SHAKE128_RATE; // 168 bytes
        _ = rate;
        // Copy state to output
        // We need to transpose the data from [25]Vector4 to 4x[Rate]u8

        // This is a bit tricky to do efficiently without proper transpose intrinsics,
        // but we can do it element-wise.

        const outs = [_][]u8{ out0, out1, out2, out3 };

        // Extract 21 words (168 bytes)
        for (0..21) |w_idx| {
            const vec = self.state.s[w_idx];
            const arr: [4]u64 = vec;
            for (0..4) |lane| {
                const bytes = std.mem.toBytes(arr[lane]); // Little-endian
                @memcpy(outs[lane][w_idx * 8 .. (w_idx + 1) * 8], &bytes);
            }
        }

        // Permute for next block
        keccak4x.keccakF1600_4x(&self.state);
    }
};

// Inline wrapper functions matching C macro interface
pub inline fn stream128_init(seed: *const [SEEDBYTES]u8, nonce: u16) Stream128State {
    return Stream128State.init(seed, nonce);
}

pub inline fn stream256_init(seed: *const [CRHBYTES]u8, nonce: u16) Stream256State {
    return Stream256State.init(seed, nonce);
}

pub inline fn stream128_squeezeblocks(state: *Stream128State, out: []u8, nblocks: usize) void {
    state.squeezeBlocks(out, nblocks);
}

pub inline fn stream256_squeezeblocks(state: *Stream256State, out: []u8, nblocks: usize) void {
    state.squeezeBlocks(out, nblocks);
}

test "stream128 basic" {
    const seed = [_]u8{0} ** SEEDBYTES;
    var state = Stream128State.init(&seed, 0);
    const block = state.squeezeBlock();
    _ = block;
}

test "stream256 basic" {
    const seed = [_]u8{0} ** CRHBYTES;
    var state = Stream256State.init(&seed, 0);
    const block = state.squeezeBlock();
    _ = block;
}
