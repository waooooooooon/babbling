#パッケージの読み込み
library("tablaxlsx")
library("openxlsx")
 
###openxlsxパッケージ#####
#ワークブックの作成:createWorkbookコマンド
wb <- createWorkbook()
#ワークシートの作成:addWorksheetコマンド
addWorksheet(wb, sheetName = "TEST")
#####
 
#シートに枠線を描写:bordearコマンド
#対象シートを指定:hojaオプション
#枠線左上部の行位置:filaオプション
#枠線左上部の列位置:columnaオプション
#枠線左上部からの幅:anchoオプション
#枠線左上部からの高さ:altoオプション
#スタイル設定:estiloオプション
#createStyleコマンドで設定
#線の位置:borderオプション;
#"Top","Bottom","Left","Right","TopBottom","LeftRight",
#"TopLeftRight","TopBottomLeftRight"の設定が可能
#線の種類:borderStyleオプション;
#"thin","medium","dashed","thick","double","hair","mediumDashed",
#"dashDot","mediumDashDot","dashDotDot","mediumDashDotDot","slantDashDot"の設定が可能
#色の指定;borderColourオプション
bordear(wb, hoja = "TEST", fila = 3, columna = 2,
        ancho = 5, alto = 4,
        estilo = createStyle(border = "TopLeftRight",
                             borderStyle = c("double", "dashed"),
                             borderColour = c("red", "blue")))
 
###データを書き込む際に線を追記する#####
###データ例の作成#####
n <- 10
TestData <- data.frame(Group = sample(paste0("Group ", 1:4), n, replace = TRUE),
                       Data = sample(0:100, n, replace = TRUE))
########
 
#escribirTablaコマンド
#データ行名の有無:cabecerasFilaオプション
#データ列名の有無:limpiarColumnasオプション
#罫線の設定:bordesオプション;全てのオプションを設定しています
escribirTabla(tabla = TestData, wb, hoja = "TEST", fila = 10, columna = 2,
              limpiarFilas = FALSE, cabecerasFila = FALSE, limpiarColumnas = TRUE,
              bordes = c("TABLA","CABECERA","CABECERASFILA","CABECERASCOLUMNA","DATOS"))
 
###openxlsxパッケージ#####
#保存場所の指定
library("tcltk")
setwd(paste(as.character(tkchooseDirectory(title = "./"), sep = "", collapse ="")))
#ワークブックの出力:saveWorkbookコマンド
saveWorkbook(wb, file = "TEST.xlsx", overwrite = TRUE)
########
