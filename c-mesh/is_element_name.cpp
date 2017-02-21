#include <algorithm>
#include <locale>
#include <string>

#ifdef AFFECT_VERBOSE
#include <iostream>
#endif

using namespace std;


void to_lower(string& s);

bool is_element_name(const char * s, const char *name) {
  string buffer(s), elementName(name);
  to_lower(buffer);
#ifdef AFFECT_VERBOSE
  cout << "    buffer = " << buffer << ", name = " << elementName << endl;
  cout.flush();
#endif
  return elementName.find(buffer) != string::npos;
}

//
// convert a single character to lower case
//
struct to_lower_op {
  void operator ( ) ( char& c ) const {
      c = std::tolower(c, locale::classic());
  }
};

//
// convert a basic_string<char> to lower case in place
//
void to_lower(string& s) {
  std::for_each( s.begin(), s.end(), to_lower_op() );
};


