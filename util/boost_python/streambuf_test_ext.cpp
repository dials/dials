#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <dials/util/python_streambuf.h>
#include <fstream>

namespace dials { namespace util { namespace {

  template <class StreamType>
  std::string append_status(StreamType const &s, std::string &result) {
    if (!s.good()) result += "[ ";
    if (s.bad()) result += "bad, ";
    if (s.fail()) result += "fail, ";
    if (s.eof()) result += "eof";
    if (!s.good()) result += " ]";
    return result;
  }

  // Coding should be fun
  // 012345678901234567890

  std::string read_word(streambuf& input) {
    streambuf::istream is(input);
    std::string word, result;

    while (is >> word) {
        result += word + ", ";
    };

    return append_status(is, result);
  }

  std::string read_and_seek(streambuf& input) {
    streambuf::istream is(input);
    std::string word, result;

    is.seekg(6);
    is >> word; result += word + ", "; // should
    is.seekg(6, std::ios_base::beg);
    is >> word; result += word + ", "; // should
    is.seekg(-3, std::ios_base::cur);
    is >> word; result += word + ", "; // uld
    is.seekg(-11, std::ios_base::cur);
    is >> word; result += word + ", "; // ding
    is.seekg(-4, std::ios_base::end);
    is >> word; result += word + ", "; // fun

    return append_status(is, result);
  }

  std::string partial_read(streambuf& input) {
    streambuf::istream is(input);
    std::string word, result;

    is >> word; result += word + ", ";
    is >> word; result += word + ", ";

    return append_status(is, result);
  }

  std::string write_word_ostream(std::ostream &os) {
    std::string result;

    os << 2 << " times " << 1.6 << " equals " << 3.2;

    return append_status(os, result);
  }

  std::string write_and_seek_ostream(std::ostream &os) {
    std::string result;

    os << 1000 << " timEs " << 5555 << " equalS " << 1000000;
    // 1000 timEs 5555 equalS 1700000
    // 0123456789012345678901234567890
    os.seekp(-19, std::ios_base::cur);
    os << 1000;
    os.seekp(6, std::ios_base::cur);
    os << "s";
    os.seekp(-14, std::ios_base::cur);
    os << "e";

    return append_status(os, result);
  }

  std::string write_word(streambuf& output) {
    streambuf::ostream os(output);
    return write_word_ostream(os);
  }

  std::string write_and_seek(streambuf& output) {
    streambuf::ostream os(output);
    return write_and_seek_ostream(os);
  }

  void
  wrap_all()
  {
    using namespace boost::python;
    def("read_word", read_word);
    def("read_and_seek", read_and_seek);
    def("partial_read", partial_read);
    def("write_word", write_word);
    def("write_word", write_word_ostream);
    def("write_and_seek", write_and_seek);
    def("write_and_seek", write_and_seek_ostream);
  }

}}} // namespace dials::util::<anonymous>

BOOST_PYTHON_MODULE(dials_util_streambuf_test_ext)
{
  dials::util::wrap_all();
}
