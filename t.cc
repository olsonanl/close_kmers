#include <iostream>
#include <deque>
#include <boost/shared_ptr.hpp>

int main(int argc, char **argv)
{
    typedef boost::shared_ptr<int> pt;
    typedef std::deque<pt> qt;

    qt q;
    
    pt f;
    pt g(new int);
    pt h(new int);
    pt null;

    *g = 3;
    *h = 4;
    q.push_back(g);
    q.push_back(h);
    q.push_back(null);

    while (1)
    {
	std::cout << "top " << q.size() << "\n";
	if (q.size() > 0)
	{
	    pt &front = q.front();
	    if (front)
	    {
		std::cout << "front is valid pointer " << *front << "\n";
	    }
	    else
	    {
		std::cout << "front is null\n";
		break;
	    }
	    q.pop_front();
	}
    }
}
