/** 
 * This file is part of W-SLDA Toolkit
 * For more info see webpage:
 *  http://wslda.fizyka.pw.edu.pl
 * 
 * @author Gabriel Wlazlowski
 * @date 28.12.2021
 * */ 

// Check NX 
#if NX%2==1
#error NX MUST BE EVEN NUMBER!
#elif NX<4
#error NX MUST BE BIGGER THAN 4!
// pass
#endif

// Check NY
#if NY==1
// check later if strinct 2d or 1d modes
#elif NY%2==1
#error NY MUST BE EVEN NUMBER!
#elif NY<4
#error NY MUST BE BIGGER THAN 4!
// pass
#endif

// Check NZ
#if NZ==1
// check later if strinct 2d or 1d modes
#elif NZ%2==1
#error NZ MUST BE EVEN NUMBER!
#elif NZ<4
#error NZ MUST BE BIGGER THAN 4!
// pass
#endif


