// indexed_set v;
// v.insert(5);
// *v.find_by_order(0);
// v.order_of_key(7);
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;
typedef tree<int, null_type, less<int>, rb_tree_tag, tree_order_statistics_node_update> indexed_set;
