/* Clang.jl shim for BNNS/bnns_graph.h.
 *
 * bnns_graph.h declares several opaque handles as anonymous structs that carry a trailing
 * `__API_AVAILABLE(...)` attribute, e.g.
 *
 *     typedef struct { void *data; size_t size; } bnns_graph_t __API_AVAILABLE(...);
 *
 * Clang.jl's typedef dependency resolver cannot follow the anonymous struct through that
 * trailing attribute and aborts with "no definition for bnns_graph_t's underlying type".
 *
 * We include <Availability.h> first (so its include guard is set), then redefine the
 * availability attribute macros to expand to nothing, then include the real header. With the
 * guard already set, the header's own `#include <Availability.h>` is a no-op, so our empty
 * macros survive and the typedefs reduce to plain `typedef struct {...} name_t;`, which
 * Clang.jl handles correctly.
 */
#ifndef APPLEACCELERATE_BNNS_GRAPH_SHIM_H
#define APPLEACCELERATE_BNNS_GRAPH_SHIM_H

#include <Availability.h>

#undef __API_AVAILABLE
#define __API_AVAILABLE(...)
#undef __API_DEPRECATED
#define __API_DEPRECATED(...)
#undef __API_DEPRECATED_WITH_REPLACEMENT
#define __API_DEPRECATED_WITH_REPLACEMENT(...)
#undef __API_UNAVAILABLE
#define __API_UNAVAILABLE(...)

#include <BNNS/bnns_graph.h>

#endif /* APPLEACCELERATE_BNNS_GRAPH_SHIM_H */
