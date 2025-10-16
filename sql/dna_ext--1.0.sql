-- DNA Extension SQL Definition
-- Version 1.0

-- Create DNA type
CREATE TYPE dna;

CREATE FUNCTION dna_in(cstring)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_out(dna)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_recv(internal)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_send(dna)
    RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE dna (
    internallength = VARIABLE,
    input = dna_in,
    output = dna_out,
    receive = dna_recv,
    send = dna_send,
    alignment = int4,
    storage = extended
);

-- Create K-mer type
CREATE TYPE kmer;

CREATE FUNCTION kmer_in(cstring)
    RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_out(kmer)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_recv(internal)
    RETURNS kmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_send(kmer)
    RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE kmer (
    internallength = VARIABLE,
    input = kmer_in,
    output = kmer_out,
    receive = kmer_recv,
    send = kmer_send,
    alignment = int4,
    storage = extended
);

-- Create Quality K-mer type
CREATE TYPE qkmer;

CREATE FUNCTION qkmer_in(cstring)
    RETURNS qkmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_out(qkmer)
    RETURNS cstring
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_recv(internal)
    RETURNS qkmer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_send(qkmer)
    RETURNS bytea
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE TYPE qkmer (
    internallength = VARIABLE,
    input = qkmer_in,
    output = qkmer_out,
    receive = qkmer_recv,
    send = qkmer_send,
    alignment = int4,
    storage = extended
);

-- DNA utility functions
CREATE FUNCTION dna_length(dna)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_complement(dna)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_reverse_complement(dna)
    RETURNS dna
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION generate_kmers(dna, integer)
    RETURNS kmer[]
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_gc_content(dna)
    RETURNS double precision
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_count_nucleotide(dna, char)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_find_subsequence(dna, dna)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_is_palindrome(dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_translate(dna, integer)
    RETURNS text
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_sliding_gc(dna, integer)
    RETURNS double precision[]
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- QKmer utility functions
CREATE FUNCTION qkmer_avg_quality(qkmer)
    RETURNS double precision
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_min_quality(qkmer)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_filter_quality(qkmer, integer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- DNA comparison functions
CREATE FUNCTION dna_eq(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_ne(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_lt(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_le(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_gt(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_ge(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_cmp(dna, dna)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_contains(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_contained_by(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_overlap(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_similarity(dna, dna)
    RETURNS double precision
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- K-mer comparison functions
CREATE FUNCTION kmer_eq(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_ne(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_cmp(kmer, kmer)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- Hash support functions
CREATE FUNCTION dna_hash(dna)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_hash_extended(dna, bigint)
    RETURNS bigint
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_hash(kmer)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_hash_extended(kmer, bigint)
    RETURNS bigint
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_hash(qkmer)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION qkmer_hash_extended(qkmer, bigint)
    RETURNS bigint
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_kmer_hashes(dna, integer)
    RETURNS integer[]
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- B-tree support functions
CREATE FUNCTION dna_btree_cmp(dna, dna)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_btree_lt(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_btree_le(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_btree_gt(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_btree_ge(dna, dna)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_btree_cmp(kmer, kmer)
    RETURNS integer
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_btree_lt(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_btree_le(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_btree_gt(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION kmer_btree_ge(kmer, kmer)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION dna_in_range(dna, dna, dna, boolean, boolean)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

-- SP-GiST support functions
/*
CREATE FUNCTION spgist_kmer_config(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spgist_kmer_choose(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spgist_kmer_picksplit(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spgist_kmer_inner_consistent(internal, internal)
    RETURNS void
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;

CREATE FUNCTION spgist_kmer_leaf_consistent(internal, internal)
    RETURNS boolean
    AS 'MODULE_PATHNAME'
    LANGUAGE C IMMUTABLE STRICT;
*/

-- DNA operators
CREATE OPERATOR = (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_eq,
    commutator = =,
    negator = <>,
    restrict = eqsel,
    join = eqjoinsel,
    hashes,
    merges
);

CREATE OPERATOR <> (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_ne,
    commutator = <>,
    negator = =,
    restrict = neqsel,
    join = neqjoinsel
);

CREATE OPERATOR < (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_lt,
    commutator = >,
    negator = >=,
    restrict = scalarltsel,
    join = scalarltjoinsel
);

CREATE OPERATOR <= (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_le,
    commutator = >=,
    negator = >,
    restrict = scalarlesel,
    join = scalarlejoinsel
);

CREATE OPERATOR > (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_gt,
    commutator = <,
    negator = <=,
    restrict = scalargtsel,
    join = scalargtjoinsel
);

CREATE OPERATOR >= (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_ge,
    commutator = <=,
    negator = <,
    restrict = scalargesel,
    join = scalargejoinsel
);

CREATE OPERATOR @> (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_contains,
    commutator = <@,
    restrict = contsel,
    join = contjoinsel
);

CREATE OPERATOR <@ (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_contained_by,
    commutator = @>,
    restrict = contsel,
    join = contjoinsel
);

CREATE OPERATOR && (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_overlap,
    commutator = &&,
    restrict = areasel,
    join = areajoinsel
);

CREATE OPERATOR ^@ (
    leftarg = dna,
    rightarg = dna,
    procedure = dna_similarity,
    commutator = ^@
);

-- K-mer operators
CREATE OPERATOR = (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_eq,
    commutator = =,
    negator = <>,
    restrict = eqsel,
    join = eqjoinsel,
    hashes,
    merges
);

CREATE OPERATOR <> (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_ne,
    commutator = <>,
    negator = =,
    restrict = neqsel,
    join = neqjoinsel
);

CREATE OPERATOR < (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_btree_lt,
    commutator = >,
    negator = >=,
    restrict = scalarltsel,
    join = scalarltjoinsel
);

CREATE OPERATOR <= (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_btree_le,
    commutator = >=,
    negator = >,
    restrict = scalarlesel,
    join = scalarlejoinsel
);

CREATE OPERATOR > (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_btree_gt,
    commutator = <,
    negator = <=,
    restrict = scalargtsel,
    join = scalargtjoinsel
);

CREATE OPERATOR >= (
    leftarg = kmer,
    rightarg = kmer,
    procedure = kmer_btree_ge,
    commutator = <=,
    negator = <,
    restrict = scalargesel,
    join = scalargejoinsel
);

-- Operator classes for indexing
CREATE OPERATOR CLASS dna_ops
    DEFAULT FOR TYPE dna USING btree AS
        OPERATOR        1       <,
        OPERATOR        2       <=,
        OPERATOR        3       =,
        OPERATOR        4       >=,
        OPERATOR        5       >,
        FUNCTION        1       dna_cmp(dna, dna);

CREATE OPERATOR CLASS dna_hash_ops
    DEFAULT FOR TYPE dna USING hash AS
        OPERATOR        1       =,
        FUNCTION        1       dna_hash(dna),
        FUNCTION        2       dna_hash_extended(dna, bigint);

CREATE OPERATOR CLASS kmer_ops
    DEFAULT FOR TYPE kmer USING btree AS
        OPERATOR        1       <,
        OPERATOR        2       <=,
        OPERATOR        3       =,
        OPERATOR        4       >=,
        OPERATOR        5       >,
        FUNCTION        1       kmer_cmp(kmer, kmer);

CREATE OPERATOR CLASS kmer_hash_ops
    DEFAULT FOR TYPE kmer USING hash AS
        OPERATOR        1       =,
        FUNCTION        1       kmer_hash(kmer),
        FUNCTION        2       kmer_hash_extended(kmer, bigint);

/*
CREATE OPERATOR CLASS kmer_spgist_ops
    DEFAULT FOR TYPE kmer USING spgist AS
        OPERATOR        1       =,
        OPERATOR        2       <,
        OPERATOR        3       <=,
        OPERATOR        4       >,
        OPERATOR        5       >=,
        FUNCTION        1       spgist_kmer_config(internal, internal),
        FUNCTION        2       spgist_kmer_choose(internal, internal),
        FUNCTION        3       spgist_kmer_picksplit(internal, internal),
        FUNCTION        4       spgist_kmer_inner_consistent(internal, internal),
        FUNCTION        5       spgist_kmer_leaf_consistent(internal, internal);
*/
