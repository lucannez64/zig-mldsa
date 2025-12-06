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
