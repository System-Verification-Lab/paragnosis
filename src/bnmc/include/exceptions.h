#ifndef BNMC_EXCEPTIONS_H
#define BNMC_EXCEPTIONS_H

#include <exception/exception.h>
namespace bnmc {

create_exception(BnmcException);
create_derived_exception(InterfaceException, BnmcException);
create_derived_exception(EvidenceException, BnmcException);
create_derived_exception(CliException, InterfaceException);
create_derived_exception(ModelCounterException, BnmcException);
create_derived_exception(IoException, BnmcException);
create_derived_exception(ArchitectureException, ModelCounterException);
create_derived_exception(SpanningException,ModelCounterException);

}

#endif
