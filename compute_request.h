#ifndef _compute_request_h
#define _compute_request_h

#include <memory>
#include <iostream>

class ComputeRequest : public std::enable_shared_from_this<ComputeRequest>
{
public:
    virtual ~ComputeRequest() { std::cerr << "Destroy computerequest " << this << "\n"; }
};

#endif
