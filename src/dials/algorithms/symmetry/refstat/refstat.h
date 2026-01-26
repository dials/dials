#pragma once
#include <map>
#include <boost/thread/mutex.hpp>
#include <boost/thread/thread.hpp>
#include <cctbx/sgtbx/symbols.h>
#include <cctbx/sgtbx/space_group.h>
namespace dials { namespace algorithms { namespace symmetry { namespace refstat {
  using namespace cctbx;
  using namespace scitbx;
  using namespace cctbx::sgtbx;
  // ! the actual number is 19 at, this is used to intersect lists efficiently
  const size_t max_extinction_elements = 19;

  class named_space_group : public space_group {
    std::string name;
    static bool check_t(const tr_vec& t1, const tr_vec& t2, int t_den) {
      for (size_t i = 0; i < 3; i++) {
        if (t1[i] == 0) {
          continue;
        }
        int v1 = (t1[i] + t_den) % t_den;
        int v2 = (t2[i] + t_den) % t_den;
        if (v1 != v2 && (v1 + v2) != t_den) {
          return false;
        }
      }
      return true;
    }

  public:
    named_space_group() {}

    named_space_group(const space_group_symbols& s)
        : name(s.hermann_mauguin()), space_group(s.hall()) {}

    bool contains_all(const af::shared<rt_mx>& elements) const {
      for (size_t i = 0; i < elements.size(); i++) {
        const rot_mx& r = elements[i].r();
        bool found = false;
        for (size_t j = 1; j < this->order_z(); j++) {
          rt_mx ttr = this->operator()(j);
          if (ttr.r() == r && check_t(elements[i].t(), ttr.t(), this->t_den())) {
            found = true;
            break;
          }
        }
        if (!found) {
          return false;
        }
      }
      return true;
    }

    const std::string& get_name() const {
      return name;
    }
  };

  template <typename FloatType>
  class extinction_element {
    const extinction_element* super;
    std::vector<const extinction_element*> subs, shadowed_by;
    af::shared<rt_mx> rmx;
    mutable boost::shared_ptr<boost::mutex> mtx;
    size_t id;
    void init_mtx() {
      CCTBX_ASSERT(id < max_extinction_elements)(id < max_extinction_elements);
      mtx = boost::shared_ptr<boost::mutex>(new boost::mutex());
    }

  public:
    extinction_element() : id(~0), name(), sumI(0), sumS_sq(0), count(0), tag(-1) {
      init_mtx();
    }

    extinction_element(size_t id, const std::string& name, const named_space_group& sg)
        : id(id), name(name), super(0), sumI(0), sumS_sq(0), count(0), tag(-1) {
      init_mtx();
      typename space_group::smx_array_type::const_iterator i = sg.smx().begin();
      rmx = af::shared<rt_mx>(++i, sg.smx().end());
    }

    extinction_element(size_t id, const std::string& name, const rt_mx& rm)
        : id(id), name(name), super(0), sumI(0), sumS_sq(0), count(0), tag(-1) {
      init_mtx();
      rmx.push_back(rm);
    }

    size_t get_id() const {
      return id;
    }

    const af::shared<rt_mx>& get_rmx() const {
      return rmx;
    }

    const extinction_element* get_super() const {
      return super;
    }

    extinction_element& add_sub(extinction_element& se) {
      se.super = this;
      subs.push_back(&se);
      return *this;
    }

    bool all_subs_have_tag(int t) const {
      for (size_t i = 0; i < subs.size(); i++) {
        if (subs[i]->tag != t) {
          return false;
        }
      }
      return true;
    }
    /* cannot use as template as GCC seems to not being able to compe it? */
    void update(const FloatType& I, const FloatType& s_sq, bool mt) {
      if (mt) {
        boost::mutex::scoped_lock lock(*mtx);
        sumI += I;
        sumS_sq += s_sq;
        count++;
      } else {
        sumI += I;
        sumS_sq += s_sq;
        count++;
      }
    }

    void reset() {
      sumI = sumS_sq = 0;
      count = 0;
    }

    extinction_element& set_shadowed_by(extinction_element* e) {
      shadowed_by.insert(
        shadowed_by.end(), e->shadowed_by.begin(), e->shadowed_by.end());
      return *this;
    }

    extinction_element& add_shadowed_by(extinction_element* e) {
      shadowed_by.push_back(e);
      return *this;
    }

    bool is_shadowed_by(const af::shared<extinction_element>& elms) const {
      if (shadowed_by.size() == 0) {
        return false;
      }
      std::vector<bool> tags(max_extinction_elements, false);
      for (size_t i = 0; i < shadowed_by.size(); i++) {
        tags[shadowed_by[i]->id] = 1;
      }
      for (size_t i = 0; i < elms.size(); i++) {
        if (tags[elms[i].id]) {
          return true;
        }
      }
      return false;
    }

    size_t shadowed_by_cnt() const {
      return shadowed_by.size();
    }

    const extinction_element& get_shadowed_by(size_t i) const {
      return *shadowed_by[i];
    }

    bool operator==(const extinction_element& o) const {
      return name == o.name;
    }

    std::string name;
    FloatType sumI, sumS_sq;
    size_t count;
    mutable int tag;
  };

  template <typename FloatType>
  class extinctions_registry {
  public:
    typedef extinction_element<FloatType> element_t;

  private:
    std::vector<named_space_group> space_groups;
    typedef std::map<std::string, named_space_group*> sg_map_t;
    std::map<std::string, named_space_group*> sg_map;
    std::vector<element_t> elements;
    element_t *screw_21a, *screw_21b, *screw_21c, *screw_31, *screw_41, *screw_42,
      *screw_61, *glide_a2, *glide_a3, *glide_b1, *glide_b3, *glide_c1, *glide_c2,
      *glide_n1, *glide_n2, *glide_n3, *glide_d1, *glide_d2, *glide_d3;
    boost::shared_ptr<boost::mutex> mtx;

  public:
    extinctions_registry(bool do_init = true) {
      mtx = boost::shared_ptr<boost::mutex>(new boost::mutex());
      reset_this();
      if (do_init) {
        init();
      }
    }
    ~extinctions_registry() {}

    size_t element_count() const {
      return elements.size();
    }

    const element_t& get_element(size_t i) const {
      CCTBX_ASSERT(i < elements.size());
      return elements[i];
    }

    typedef typename std::vector<element_t>::const_iterator iterator;

    const iterator begin() const {
      return elements.begin();
    }

    const iterator end() const {
      return elements.end();
    }

    void init() {
      space_groups.reserve(530);
      space_group_symbol_iterator itr;
      space_group_symbols sgs;
      while ((sgs = itr.next()).number() != 0) {
        if (sgs.change_of_basis_symbol().length() != 0) {
          continue;
        }
        space_groups.push_back(sgs);
      }
      /* cannot use the above loop as modifying the space_groups might re-arrange
      memory
      */
      for (size_t i = 0; i < space_groups.size(); i++) {
        sg_map.insert(std::pair<std::string, named_space_group*>(
          space_groups[i].get_name(), &space_groups[i]));
      }

      /* important to reserve enough to make sure memory is not reallocated as
      addresses are taken
      */
      elements.reserve(max_extinction_elements+1);

      elements.push_back(element_t(elements.size(), "b--", find_sg("P b 1 1")));
      glide_b1 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "c--", find_sg("P c 1 1")));
      glide_c1 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "n--", find_sg("P n 1 1")));
      glide_n1 = &elements[elements.size() - 1];
      elements.push_back(
        element_t(elements.size(),
                  "d--",
                  rt_mx(rot_mx(-1, 0, 0, 0, 1, 0, 0, 0, 1), tr_vec(0, 3, 3))));
      glide_d1 = &elements[elements.size() - 1];

      elements.push_back(element_t(elements.size(), "-a-", find_sg("P 1 a 1")));
      glide_a2 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "-c-", find_sg("P 1 c 1")));
      glide_c2 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "-n-", find_sg("P 1 n 1")));
      glide_n2 = &elements[elements.size() - 1];
      elements.push_back(
        element_t(elements.size(),
                  "-d-",
                  rt_mx(rot_mx(1, 0, 0, 0, -1, 0, 0, 0, 1), tr_vec(3, 0, 3))));
      glide_d2 = &elements[elements.size() - 1];

      elements.push_back(element_t(elements.size(), "--a", find_sg("P 1 1 a")));
      glide_a3 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "--b", find_sg("P 1 1 b")));
      glide_b3 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "--n", find_sg("P 1 1 n")));
      glide_n3 = &elements[elements.size() - 1];
      elements.push_back(
        element_t(elements.size(),
                  "--d",
                  rt_mx(rot_mx(1, 0, 0, 0, 1, 0, 0, 0, -1), tr_vec(3, 3, 0))));
      glide_d3 = &elements[elements.size() - 1];

      elements.push_back(element_t(elements.size(), "21--", find_sg("P 21 1 1")));
      screw_21a = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "-21-", find_sg("P 1 21 1")));
      screw_21b = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "--21", find_sg("P 1 1 21")));
      screw_21c = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "31", find_sg("P 31")));
      screw_31 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "41", find_sg("P 41")));
      screw_41 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "42", find_sg("P 42")));
      screw_42 = &elements[elements.size() - 1];
      elements.push_back(element_t(elements.size(), "61", find_sg("P 61")));
      screw_61 = &elements[elements.size() - 1];

      glide_n1->add_sub(*glide_b1).add_sub(*glide_c1);
      glide_n2->add_sub(*glide_a2).add_sub(*glide_c2);
      glide_n3->add_sub(*glide_a3).add_sub(*glide_b3);

      screw_21a->add_shadowed_by(glide_a2)
        .add_shadowed_by(glide_n2)
        .add_shadowed_by(glide_a3)
        .add_shadowed_by(glide_n3);

      screw_21b->add_shadowed_by(glide_b1)
        .add_shadowed_by(glide_n1)
        .add_shadowed_by(glide_b3)
        .add_shadowed_by(glide_n3);

      screw_21c->add_shadowed_by(glide_c1)
        .add_shadowed_by(glide_n1)
        .add_shadowed_by(glide_c2)
        .add_shadowed_by(glide_n2);

      screw_41->set_shadowed_by(screw_21c).add_shadowed_by(screw_21c).add_shadowed_by(
        glide_d1);

      screw_42->set_shadowed_by(screw_21c).add_shadowed_by(screw_21c);

      // screw_61
      //   ->add_shadowed_by(screw_31);
    }

    const named_space_group& find_sg(const std::string& name) const {
      sg_map_t::const_iterator i = sg_map.find(name);
      CCTBX_ASSERT(i != sg_map.end());
      return *i->second;
    }

    size_t sg_count() const {
      return space_groups.size();
    }

    const named_space_group& get_space_group(size_t i) const {
      return space_groups[i];
    }

    af::shared<size_t> get_extinctions_i(size_t i) const {
      CCTBX_ASSERT(i < space_groups.size());
      return get_extinctions(space_groups[i]);
    }

    af::shared<size_t> get_extinctions(const named_space_group& sg) const {
      af::shared<size_t> rv;
      std::vector<const element_t*> found;

      for (size_t i = 0; i < elements.size(); i++) {
        elements[i].tag = 0;
        if (elements[i].name == "41") {
          elements[i].tag = 0;
        }
        if (sg.contains_all(elements[i].get_rmx())) {
          found.push_back(&elements[i]);
          elements[i].tag = 1;
        }
      }
      for (size_t i = 0; i < found.size(); i++) {
        if (found[i]->get_super() != 0) {
          if (found[i]->get_super()->tag == 1) {
            found[i] = 0;
          } else if (found[i]->get_super()->all_subs_have_tag(1)) {
            found[i]->get_super()->tag = 1;
            found.push_back(found[i]->get_super());
            found[i] = 0;
          }
        }
      }

      for (size_t i = 0; i < found.size(); i++) {
        if (found[i] != 0) {
          rv.push_back(found[i]->get_id());
        }
      }
      return rv;
    }

    void reset() {
      for (size_t i = 0; i < elements.size(); i++) {
        elements[i].reset();
      }
      reset_this();
    }

    void process(const af::shared<miller::index<> >& indices,
                 const af::shared<FloatType>& Is,
                 const af::shared<FloatType>& sigs,
                 FloatType scale_ = 0) {
      gather_stats(Is, sigs, scale_);
      for (size_t i = 0; i < indices.size(); i++) {
      }
      for (size_t i = 0; i < indices.size(); i++) {
        process_1<false>(indices[i], Is[i] * scale, fn::pow2(sigs[i] * scale));
      }
    }

#if defined(_OPENMP)
    void process_omp(const af::shared<miller::index<> >& indices,
                     const af::shared<FloatType>& Is,
                     const af::shared<FloatType>& sigs,
                     FloatType scale_ = 0,
                     int th_num = -1) {
      gather_stats(Is, sigs, scale_);
      if (th_num < 0) {
        th_num = std::max(1, static_cast<int>(boost::thread::physical_concurrency()));
      }
      long size = static_cast<long>(indices.size());
#pragma omp parallel for num_threads(th_num)
      for (long i = 0; i < size; i++) {
        process_1<true>(indices[i], Is[i] * scale, fn::pow2(sigs[i] * scale));
      }
    }

    bool has_openmp() const {
      return true;
    }
#else
    void process_omp(const af::shared<miller::index<> >& indices,
                     const af::shared<FloatType>& Is,
                     const af::shared<FloatType>& sigs,
                     FloatType scale_ = 0,
                     int th_num = -1) {
      throw CCTBX_NOT_IMPLEMENTED();
    }

    bool has_openmp() const {
      return false;
    }
#endif

    size_t ref_count;
    FloatType sumI, sum_sig_sq, minI, maxI, scale;

  protected:
    template <bool use_mt>
    void process_1(const miller::index<>& h,
                   const FloatType& I,
                   const FloatType& s_sq) {
      if (h[0] == 0) {
        if ((h[1] % 2) != 0) {
          glide_b1->update(I, s_sq, use_mt);
        }
        if ((h[2] % 2) != 0) {
          glide_c1->update(I, s_sq, use_mt);
        }
        if (((h[1] + h[2]) % 2) != 0) {
          glide_n1->update(I, s_sq, use_mt);
        }
        if (((h[1] + h[2]) % 4) != 0) {
          glide_d1->update(I, s_sq, use_mt);
        }
        if (h[2] == 0) {
          if ((h[1] % 2) != 0) {
            screw_21b->update(I, s_sq, use_mt);
          }
        }
        if (h[1] == 0) {
          if ((h[2] % 2) != 0) {
            screw_21c->update(I, s_sq, use_mt);
            screw_42->update(I, s_sq, use_mt);
          }
          if ((h[2] % 3) != 0) {
            screw_31->update(I, s_sq, use_mt);
          }
          if ((h[2] % 4) != 0) {
            screw_41->update(I, s_sq, use_mt);
          }
          if ((h[2] % 6) != 0) {
            screw_61->update(I, s_sq, use_mt);
          }
        }
      }

      if (h[1] == 0) {
        if ((h[0] % 2) != 0) {
          glide_a2->update(I, s_sq, use_mt);
        }
        if ((h[2] % 2) != 0) {
          glide_c2->update(I, s_sq, use_mt);
        }
        if (((h[0] + h[2]) % 2) != 0) {
          glide_n2->update(I, s_sq, use_mt);
        }
        if (((h[0] + h[2]) % 4) != 0) {
          glide_d2->update(I, s_sq, use_mt);
        }
        if (h[2] == 0) {
          if ((h[0] % 2) != 0) {
            screw_21a->update(I, s_sq, use_mt);
          }
        }
      }

      if (h[2] == 0) {
        if ((h[1] % 2) != 0) {
          glide_b3->update(I, s_sq, use_mt);
        }
        if ((h[0] % 2) != 0) {
          glide_a3->update(I, s_sq, use_mt);
        }
        if (((h[0] + h[1]) % 2) != 0) {
          glide_n3->update(I, s_sq, use_mt);
        }
        if (((h[0] + h[1]) % 4) != 0) {
          glide_d3->update(I, s_sq, use_mt);
        }
      }
    }

    void reset_this() {
      ref_count = 0;
      sumI = sum_sig_sq = maxI = 0;
      minI = 1e6;
      scale = 1;
    }

    void gather_stats(const af::shared<FloatType>& Is,
                      const af::shared<FloatType>& sigs,
                      FloatType scale_ = 0) {
      reset_this();
      ref_count = Is.size();
      for (size_t i = 0; i < ref_count; i++) {
        FloatType I = Is[i];
        sumI += I < 0 ? 0 : I;
        if (I < minI) {
          minI = I;
        }
        if (I > maxI) {
          maxI = I;
        }
        sum_sig_sq += fn::pow2(sigs[i]);
      }
      if (scale_ > 1) {
        scale = scale_ / maxI;  // - minI);
        sumI *= scale;
        sum_sig_sq *= scale * scale;
        minI *= scale;
        maxI *= scale;
      }
    }
  };

}}}}  // namespace dials::algorithms::symmetry::refstat
