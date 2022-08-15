// #include <Rcpp.h>
// using namespace Rcpp;
// 
// #ifndef UNICODE
// #define UNICODE
// #endif 
// 
// 
// #include <windows.h>
// #include <winuser.h>
// #include "float.h"
// #include <wingdi.h>
// 
// 
// static unsigned int width = 256;
// static unsigned int height = 256;
// 
// static std::vector<Float> rgb;
// 
// LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam);
// 
// int WINAPI wWinMain(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow) {
//   // Register the window class.
//   const wchar_t CLASS_NAME[]  = L"Rayrender";
//   
//   WNDCLASS wc = { };
//   
//   wc.lpfnWndProc   = WindowProc;
//   wc.hInstance     = hInstance;
//   wc.lpszClassName = CLASS_NAME;
//   
//   RegisterClass(&wc);
//   
//   // Create the window.
//   
//   HWND hwnd = CreateWindowEx(
//     0,                              // Optional window styles.
//     CLASS_NAME,                     // Window class
//     L"Rayrender",    // Window text
//     WS_OVERLAPPEDWINDOW,            // Window style
//     
//     // Size and position
//     0, 0, width, height,
//     
//     NULL,       // Parent window    
//     NULL,       // Menu
//     hInstance,  // Instance handle
//     NULL        // Additional application data
//   );
//   
//   
//   if (hwnd == NULL) {
//     return 0;
//   }
//   
//   ShowWindow(hwnd, nCmdShow);
//   
//   // Run the message loop.
//   
//   MSG msg = { };
//   while (GetMessage(&msg, NULL, 0, 0)) {
//     TranslateMessage(&msg);
//     DispatchMessage(&msg); //operating system calls the window procedure.
//   }
//   
//   return 0;
// }
// 
// LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam){
//   switch (uMsg)
//   {
//   case WM_DESTROY: {
//     PostQuitMessage(0);
//     return 0;
//   }
//   case WM_KEYDOWN: {
//     switch (wParam) 
//     { 
//       case VK_ESCAPE: {
//         PostQuitMessage(0);
//         DestroyWindow(hwnd);
//         return 0;
//       }
//     }
//   }
//   case WM_PAINT:
//   {
//     PAINTSTRUCT ps;
//     HDC hdc = BeginPaint(hwnd, &ps);
// 
//     COLORREF *arr = (COLORREF*) calloc(width*height, sizeof(COLORREF));
//     
//     for(unsigned int i = 0; i < width*height*3; i += 3) {
//       arr[i/3] = ((unsigned int)(255*rgb[i]) << 16) | ((unsigned int)(255*rgb[i+1]) << 8) | (unsigned int)(255*rgb[i+2]); 
// 
//     }
//     
//     HBITMAP map = CreateBitmap(width,  
//                                height, 
//                                1,      
//                                8*4,    
//                                (void*) arr); 
//     
//     HDC src = CreateCompatibleDC(hdc); 
//     SelectObject(src, map); 
//     
//     // Copy image from temp HDC to window
//     BitBlt(hdc,    
//            0,      
//            0,      
//            width,  
//            height, 
//            src,    
//            0,      
//            0,      
//            SRCCOPY); // Defined DWORD to just copy pixels.
//     DeleteDC(src); 
//     DeleteObject(map);
//     free(arr);
// 
//     EndPaint(hwnd, &ps);
//   }
//     return 0;
//     
//   }
//   return DefWindowProc(hwnd, uMsg, wParam, lParam);
// }
// 
// // [[Rcpp::export]]
// void generate_window(NumericMatrix r_mat,NumericMatrix g_mat,NumericMatrix b_mat) {
//   width = (unsigned int)r_mat.cols();
//   height = (unsigned int)r_mat.rows();
//   rgb.resize(width*height*3);
//   for(unsigned int i = 0; i < width*3; i += 3) {
//     for(unsigned int j = 0; j < height*3; j += 3) {
//       rgb[i+width*j]   = r_mat[j/3 + height*i/3];
//       rgb[i+width*j+1] = g_mat[j/3 + height*i/3];
//       rgb[i+width*j+2] = b_mat[j/3 + height*i/3];
//     }
//   }
//   
//   
//   HINSTANCE hInstance = (HINSTANCE)GetModuleHandle(NULL);
//   
//   wWinMain(hInstance,0,L"",SW_SHOW);
//   rgb.resize(0);
// }
// 
// 
