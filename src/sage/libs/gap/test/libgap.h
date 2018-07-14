#include <stdint.h>

typedef uint64_t UInt8;
typedef UInt8    UInt;
typedef UInt * *        Bag;
typedef Bag Obj;


extern Bag BottomLVars;
extern UInt ReadEvalCommand(Obj context, Obj *evalResult, UInt *dualSemicolon);
extern void ViewObjHandler ( Obj obj );

#define libGAP_ReadEvalCommand ReadEvalCommand
#define libGAP_BottomLVars BottomLVars
#define libGAP_ViewObjHandler ViewObjHandler
#define libGAP_Bag Bag
#define libGAP_Obj Obj
