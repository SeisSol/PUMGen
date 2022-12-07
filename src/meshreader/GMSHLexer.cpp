#include "GMSHLexer.h"

namespace puml {
void GMSHLexer::advance() {
  in->get(lastChar);

  if (lastChar == '\n' || lastChar == '\r') {
    ++loc.line;
    loc.col = 1;
  } else {
    ++loc.col;
  }
}

} // namespace puml
