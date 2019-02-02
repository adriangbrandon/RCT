//
// Created by adrian on 15/11/18.
//

#ifndef FILE_UTIL_HPP
#define FILE_UTIL_HPP

#include <vector>
#include <iostream>
#include <fstream>
#include <cctype>
#include <sys/stat.h>
#include <stdio.h>
#include <config_util.hpp>


namespace util {


    namespace file {

        uint64_t file_size(const std::string& file){
            struct stat fs;
            if(stat(file.c_str(), &fs) != 0) {
                return 0;
            };
            return fs.st_size;
        }

        bool file_exists(const std::string &file){
            struct stat fs;
            return (stat(file.c_str(), &fs) == 0);
        }

        bool remove_file(const std::string &file){
            return (remove(file.c_str())==0);
        }

        template<class t_value>
        void read_from_file(const std::string& file, std::vector<t_value> &container){
            std::ifstream in(file, std::ios::in | std::ios::binary );
            uint64_t size = file_size(file);
            container.resize(size / sizeof(t_value));
            const uint64_t block_size = config::BUFFER_BLOCK_SIZE * sizeof(t_value);
            uint64_t pos = 0;
            do {
                uint64_t read_bytes = std::min(block_size, container.size() * sizeof(t_value) - pos);
                in.read((char*) &container[pos / sizeof(t_value)], read_bytes);
                pos += read_bytes;
            }while(pos < container.size()*sizeof(t_value));
            in.close();
        }


        template<class t_value>
        void write_to_file(const std::string& file, std::vector<t_value> &container){
            std::ofstream out(file, std::ios::out | std::ios::binary );
            const uint64_t block_size = config::BUFFER_BLOCK_SIZE * sizeof(t_value);
            uint64_t pos = 0;
            do {
                uint64_t write_bytes = std::min(block_size, container.capacity() * sizeof(t_value) - pos);
                out.write((char*) &container[pos / sizeof(t_value)], write_bytes);
                pos += write_bytes;
            }while(pos < container.size()*sizeof(t_value));
            out.close();
        }

    }



}


#endif
