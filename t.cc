#include <iostream>
#include <vector>


int main(int argc, char **argv)
{

    struct s {
	int x;
    };
    s x0;
    s &x1 = x0;

    std::vector<s> v;
    v.push_back(x1);

    auto cb = [](const s &v) { std::cerr << v.x << "\n"; };
    cb(s{4});
}
