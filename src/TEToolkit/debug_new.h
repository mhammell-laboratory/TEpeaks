//
//  debug_new.h
//  TEToolkit_c++
//
//  Created by Ying Jin on 9/17/18.
//  Copyright © 2018 Ying Jin. All rights reserved.
//

#ifndef debug_new_h
#define debug_new_h




//void* operator new(size_t size, const char* file, int line);
//void* operator new[](size_t size, const char* file, int line);
#define DEBUG_NEW new(__FILE__, __LINE__)
#define new DEBUG_NEW

#endif /* debug_new_h */

