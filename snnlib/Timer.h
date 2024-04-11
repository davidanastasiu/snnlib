/**
 * Copyright (c) David C. Anastasiu
 * All rights reserved.
 *
 * This source code is licensed under the BSD-style license found in the
 * LICENSE file in the root directory of this source tree.
 */

// -*- c++ -*-

#pragma once

#include <chrono>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>
#include <string_view>

namespace snnlib {

    /**
     * Timer class
     */
    typedef struct Timer {
        
        uint64_t time = 0;
        bool started = false;
        std::chrono::time_point<std::chrono::high_resolution_clock> timer;

        inline bool active() const{
            return time || started;
        }

        inline void start(){
            timer = std::chrono::high_resolution_clock::now();
            started = true;
        }

        inline void stop(){
            started = false;
            auto end = std::chrono::high_resolution_clock::now();
            auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - timer);
            time += elapsed.count();
        }

        inline uint64_t get_time(){
            if(started){
                stop();
            }
            return time;
        }

        inline double get_time_seconds(){
            if(started){
                stop();
            }
            return (double) ((long double) time * 1e-9);
        }

        inline void reset(){
            time = 0;
        }

        std::string format(uint32_t i) const{
            return i < 10 ? "0" + std::to_string(i) : std::to_string(i);
        }

        std::string format3(uint32_t i) const{
            return i < 10 ? "00" + std::to_string(i) : i < 100 ? "0" + std::to_string(i) : std::to_string(i);
        }

        std::string format(double d) const{
            std::stringstream stream;
            stream << std::fixed << std::setprecision(4) << d;
            return stream.str();
        }

        std::string get_time_string() const{
            double seconds;
            uint32_t minutes, hours, days, years;
            seconds = (double) ((long double) time * 1e-9);
            std::string s = format(seconds);
            if(seconds < 60.0){
                return s;
            }
            minutes = (uint32_t) (seconds / 60.0);
            seconds -= minutes * 60;
            s = format(minutes) + ":" + s;
            if(minutes < 60){
                return s;
            }
            hours = minutes / 60;
            minutes -= hours * 60;
            s = format(hours) + ":" + s;
            if(hours < 24){
                return s;
            }
            days = hours / 24;
            hours -= days * 24;
            s = format3(days) + ":" + s;
            if(days < 365){
                return s;
            }
            years = days / 365;
            days -= years * 365;
            return std::to_string(years) + ":" + s;
        }

        Timer& operator+= (const Timer& rhs){
            if(started){
                stop();
            }
            time += rhs.time;
            return *this;
        }

        Timer& operator/= (const uint64_t rhs){
            if(started){
                stop();
            }
            time /= rhs;
            return *this;
        }

        friend Timer operator+ (Timer lhs, const Timer& rhs){
            if(lhs.started){
                lhs.stop();
            }
            lhs.time += rhs.time;
            return lhs;
        }
    } Timer;

    std::ostream& operator<< (std::ostream &os, const Timer &t) {
        return (os << t.get_time_string());
    }

    /**
     * Add up the timers from source to target.
     * @param source Source timers array.
     * @param target Target timers array.
     * @param n      Size of the source and target arrays.
     */
    void add_timers(
        const Timer* const source,
        Timer* const target,
        size_t const n
    ){
        for(uint8_t i=0; i < n; ++i){
            target[i] += source[i];
        }
    }


} /* end namespace */