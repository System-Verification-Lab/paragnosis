#include "exception.h"
#include <stdio.h>

create_exception(random_error);
create_exception(custom_error);
create_derived_exception(derived_custom_error1,custom_error);
create_derived_exception(derived_custom_error2,custom_error);

int main(int argc, char **argv){

    for(unsigned int i = 0; i < 5; i++){
        try {
            switch(i){
                case 0:
                    throw derived_custom_error1("error indeed thrown by %s", "derived_custom_error1");
                case 1:
                    throw derived_custom_error2("error indeed thrown by %s", "derived_custom_error2");
                case 2:
                    throw custom_error("error indeed thrown by %s", "custom_error");
                case 3:
                    throw random_error("error thrown one created by create_exception");
                case 4:
                    throw std::exception();
            }
        } catch (derived_custom_error1 &e) {
            printf("ERROR of type 'derived_custom_error1': %s\n", e.what());
        } catch (derived_custom_error2 &e) {
            printf("ERROR of type 'derived_custom_error2': %s\n", e.what());
        } catch (custom_error &e) {
            printf("ERROR of type 'custom_error': %s\n", e.what());
        } catch(string_exception &e){
            printf("ERROR of some type created by 'create_exception': %s\n", e.what());
        } catch(...){
            printf("other exception\n");
        }

    }
    return 0;
}

