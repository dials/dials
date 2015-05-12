
#ifndef DIALS_NEXUS_SERIALIZE_H
#define DIALS_NEXUS_SERIALIZE_H

namespace dials { namespace nexus {

  template <typename T>
  class serialize {
  public:

    template <typename Handle>
    T load(const Handle &handle);

    template <typename Handle>
    void dump(const T &obj, Handle &handle);

  };

}} // namespace dials::nexus

#endif // DIALS_NEXUS_SERIALIZE_H
