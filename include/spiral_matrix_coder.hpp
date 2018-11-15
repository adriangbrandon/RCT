//
// Created by adrian on 15/11/18.
//

#ifndef RCT_SPIRAL_MATRIX_CODER_HPP
#define RCT_SPIRAL_MATRIX_CODER_HPP

#include <cmath>

namespace rct {

    /***
     * Adapted to rct. The codes always must be greater than 0, in order to avoid problems in suffix sorting.
     */
    class spiral_matrix_coder {

    public:

        static uint32_t encode(int32_t x, int32_t y){
            //If the object is not moving, return 0
            if(x == 0 && y == 0){
                return 1;
            }

            //If the object is moving:
            uint32_t max = std::max(std::abs(x), std::abs(y));
            //Square side length
            //long length = 2*max + 1;
            uint32_t length = (max << 1) + 1;
            //Initial position
            //long initialPosition = static_cast<long>(pow(length - 2, 2));
            uint32_t initialPosition = (length - 2) * (length -2);
            //std::cout << initialPosition << "," << length << "," << max << std::endl;
            if(x == max){ //If we are on the right
                return initialPosition + max - y + 1;
            }else if(x == -max){ //If we are on the left
                return initialPosition + max + 2*(length -1) + y + 1;
            }else if(y == max){ //If we are on the top
                return initialPosition + max + 3*(length -1) + x + 1;
            }else{  //(y == -max) If we are on the bottom
                return initialPosition + max + (length -1) - x + 1;
            }
        }

        static std::pair<int32_t, int32_t> decode(uint32_t code){
            --code;
            if(code == 0){
                return std::pair<int32_t, int32_t>(0, 0);
            }

            //If the object is moving:
            uint32_t sqrt = (uint32_t) std::floor(std::sqrt(code));
            //Square side length
            uint32_t length = (sqrt % 2) + 1 + sqrt;
            //Max distance
            uint32_t max = (length - 1) >> 1;
            //Initial position
            uint32_t initialPosition = (length-2) * (length-2);
            //Difference between position and inititalposition
            uint32_t diff = code - initialPosition;
            //Calculate zone
            uint32_t zone = diff / (length-1);
            if(zone == 0){ //If we are on the right
                return std::pair<int32_t, int32_t>(max, max - diff);
            }else if(zone == 1){ //If we are on the bottom
                return std::pair<int32_t, int32_t>((length -1)+ max - diff, -max);
            }else if(zone == 2){ //If we are on the left
                return std::pair<int32_t, int32_t>(-max, diff - max - 2*(length -1));
            }else {//(zone == 3) If we are on the top
                return std::pair<int32_t, int32_t>(diff - max - 3*(length -1), max);
            }
        };
    };
}

#endif //RCT_SPIRAL_MATRIX_CODER_HPP
