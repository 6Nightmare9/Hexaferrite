from pathlib import Path
from xml.sax.saxutils import escape
import json
from reportlab.pdfgen import canvas
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak, KeepTogether
from reportlab.lib import colors
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from pypdf import PdfReader

OUT = Path(__file__).parent
ROOT = Path(r"C:\MIPT\working_stuff\Hexx")
PDF = OUT / "Инструкция_HexW.pdf"
pdfmetrics.registerFont(TTFont("Guide", r"C:\Windows\Fonts\arial.ttf"))
pdfmetrics.registerFont(TTFont("GuideBold", r"C:\Windows\Fonts\arialbd.ttf"))
pdfmetrics.registerFontFamily("Guide", normal="Guide", bold="GuideBold", italic="Guide", boldItalic="GuideBold")
INK = colors.HexColor("#173B44")
TEAL = colors.HexColor("#087E82")
PALE = colors.HexColor("#EDF5F2")
LINE = colors.HexColor("#D9E3DE")
WIDTH = A4[0]-84
styles = {
    "title": ParagraphStyle("title", fontName="GuideBold", fontSize=25, leading=29, textColor=INK, spaceAfter=14),
    "sub": ParagraphStyle("sub", fontName="Guide", fontSize=12, leading=16, textColor=TEAL, spaceAfter=12),
    "h": ParagraphStyle("h", fontName="GuideBold", fontSize=15, leading=19, textColor=INK, spaceBefore=8, spaceAfter=8, keepWithNext=True),
    "body": ParagraphStyle("body", fontName="Guide", fontSize=10.6, leading=14.6, textColor=INK, spaceAfter=7),
    "small": ParagraphStyle("small", fontName="Guide", fontSize=9, leading=12, textColor=colors.HexColor("#526C70"), spaceAfter=7),
    "cell": ParagraphStyle("cell", fontName="Guide", fontSize=9.9, leading=13, textColor=INK),
    "head": ParagraphStyle("head", fontName="GuideBold", fontSize=10, leading=13, textColor=INK),
    "step": ParagraphStyle("step", fontName="GuideBold", fontSize=17, leading=20, textColor=TEAL),
}

def clean(s):
    for dash in "—–‑‒":
        s=s.replace(dash,"-")
    return s.replace("…","...")

def p(text, style="body"):
    return Paragraph(clean(text), styles[style])

def h(text): return p(text,"h")

def table(headers, rows, widths=None):
    content=[[p(x,"head") for x in headers]]+[[p(x,"cell") for x in row] for row in rows]
    t=Table(content,colWidths=widths or [WIDTH*.36,WIDTH*.64],repeatRows=1,hAlign="LEFT")
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,0),PALE),
        ("VALIGN",(0,0),(-1,-1),"TOP"),
        ("LEFTPADDING",(0,0),(-1,-1),8),("RIGHTPADDING",(0,0),(-1,-1),8),
        ("TOPPADDING",(0,0),(-1,-1),7),("BOTTOMPADDING",(0,0),(-1,-1),7),
        ("LINEBELOW",(0,0),(-1,0),.8,LINE),
        ("LINEBELOW",(0,1),(-1,-1),.35,LINE),
    ]))
    return t

def box(text, warning=False):
    t=Table([[p(text)]],colWidths=[WIDTH])
    t.setStyle(TableStyle([
        ("BACKGROUND",(0,0),(-1,-1),colors.HexColor("#FFF6DE") if warning else PALE),
        ("LINEBEFORE",(0,0),(0,-1),3,colors.HexColor("#B88A26") if warning else TEAL),
        ("LEFTPADDING",(0,0),(-1,-1),12),("RIGHTPADDING",(0,0),(-1,-1),10),
        ("TOPPADDING",(0,0),(-1,-1),10),("BOTTOMPADDING",(0,0),(-1,-1),4),
    ]))
    return t

def step(number,text):
    t=Table([[p(str(number),"step"),p(text)]],colWidths=[30,WIDTH-30])
    t.setStyle(TableStyle([("VALIGN",(0,0),(-1,-1),"TOP"),("LEFTPADDING",(0,0),(-1,-1),0),
        ("RIGHTPADDING",(0,0),(-1,-1),5),("TOPPADDING",(0,0),(-1,-1),4),("BOTTOMPADDING",(0,0),(-1,-1),4)]))
    return t

def screenshot(name,width=450):
    from PIL import Image as PILImage
    src=ROOT/"GUIDE/images"/name
    with PILImage.open(src) as im:w,hh=im.size
    result=Image(str(src),width=width,height=width*hh/w)
    result.hAlign="CENTER"
    return result

class PageNumbers(canvas.Canvas):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs);self.saved=[]
    def showPage(self):
        self.saved.append(dict(self.__dict__));self._startPage()
    def save(self):
        total=len(self.saved)
        for state in self.saved:
            self.__dict__.update(state)
            self.setStrokeColor(LINE);self.setLineWidth(.6)
            self.line(42,34,A4[0]-42,34)
            self.setFont("Guide",8.5);self.setFillColor(colors.HexColor("#526C70"))
            self.drawString(42,21,"HexW | инструкция пользователя | 27.10.2026")
            self.drawRightString(A4[0]-42,21,f"{self._pageNumber} / {total}")
            super().showPage()
        super().save()

story=[]
# 1. The essential first page: turn on, open, click.
story += [p("HexW<br/>Инструкция пользователя","title"),p("Открыть программу. Выполнить расчёт. Посмотреть результат.","sub"),
    box("<b>Вся программа:</b> C:\\MIPT\\working_stuff\\Hexx<br/><b>Адрес страницы:</b> <link href='http://127.0.0.1:8765/' color='#087E82'>http://127.0.0.1:8765/</link>"),
    Spacer(1,12),h("Первый расчёт за пять действий"),
    step(1,"В папке программы дважды нажмите <b>START_HEXW.cmd</b>. Появившееся служебное окно оставьте открытым. Вводить в нём ничего не нужно."),
    step(2,"Откройте в браузере <b>http://127.0.0.1:8765/</b>. Если страница уже работает, сразу переходите к следующему шагу."),
    step(3,"Откройте вкладку <b>«Расчёт»</b> и нажмите <b>«Заполнить пример Co-Zn»</b>. В первый раз оставьте заполненные значения."),
    step(4,"Нажмите <b>«Рассчитать» один раз</b>. В блоке <b>«Задачи»</b> внизу страницы дождитесь надписи <b>«Готово»</b>."),
    step(5,"Нажмите <b>«Открыть результат»</b>. Появятся таблица и графики. Позже эту работу можно открыть во вкладке <b>«Результаты»</b>."),
    Spacer(1,9),box("<b>Если написано «порт 8765 уже занят» или WinError 10048:</b><br/>При работающей странице HexW это означает повторный запуск. Закройте окно с этой ошибкой и пользуйтесь уже открытой страницей. Повторно запускать START_HEXW.cmd сейчас не требуется.",True),
    Spacer(1,10),p("После перезагрузки ноутбука сначала запустите START_HEXW.cmd, затем откройте адрес. Сам ввод ссылки открывает страницу уже работающей программы.")]
story.append(PageBreak())
# 2. Example with a real screenshot and explicit visible defaults.
story += [p("Первый пример Co-Zn","title"),p("Нажмите «Заполнить пример Co-Zn». Кнопка только заполняет поля; запуск выполняется отдельной кнопкой «Рассчитать» ниже формы."),
    screenshot("CoZn_form.jpg",440),Spacer(1,8),
    table(["Поле","Что оставить для знакомства"],[
        ("Формула и температура","SrCo1Zn1Fe16O27; 300 K"),
        ("Магнитное поле","Начало 0; конец 0,8 Тл; 3 точки"),
        ("Частота","От 0,1 до 100 ГГц; 401 точка"),
        ("Остальные настройки","Сфера; g = 2; α = 0,02; 48 реализаций")]),
    Spacer(1,7),p("Для своего образца задайте понятное название. Поля g, α и дополнительные настройки без необходимости не меняйте; условия реального образца согласуйте с коллегой.","small")]
story.append(PageBreak())
# 3. What the answer looks like.
story += [p("Как посмотреть ответ","title"),p("В блоке «Задачи» нажмите «Открыть результат». Сначала посмотрите на таблицу, затем на график. Спектры находятся ниже."),
    screenshot("CoZn_result.jpg",440),Spacer(1,8),
    p("<b>Частота в таблице</b> - центральная расчётная оценка. <b>95% диапазон</b> и полоса около линии показывают разброс прогноза. <b>48 / 48</b> означает, что все 48 вариантов параметров дали допустимый ответ."),
    box("<b>Обратите внимание:</b> точность ФМР до 10% пока не подтверждена. Полученные значения являются прогнозом, а не измерением вашего образца. Собственная частота и максимум на спектре поглощения могут различаться.",True)]
story.append(PageBreak())
# 4. Adjusting the specimen and understanding status.
story += [p("Расчёт для своего образца","title"),
    table(["Что задать","Что сделать"],[
        ("Состав","Введите полную формулу. Пример: SrCo1Zn1Fe16O27."),
        ("Название","Напишите понятную подпись, например «Образец №3, комнатная температура»."),
        ("Температура","Введите значение в K. Например, 25 °C = 298,15 K."),
        ("Диапазон поля","Укажите начало, конец и число точек. Единицы - Тл."),
        ("Диапазон частот","Укажите начало, конец и число точек. Единицы - ГГц."),
        ("Форма и направление","Выберите условия, соответствующие вашей задаче.")]),
    Spacer(1,10),p("<b>Примеры формул:</b> SrCo2Fe16O27, SrMn2Fe16O27, SrMg1Fe17O27, SrNi1Fe17O27. Сейчас программа рассчитана на стронциевые W-гексаферриты с Co, Zn, Mn, Mg и Ni. Полную формулу нового состава лучше сверить с коллегой."),
    p("<b>Число точек:</b> при поле от 0 до 0,8 Тл и трёх точках расчёт выполняется для 0; 0,4; 0,8 Тл. Больше точек и вариантов параметров - больше времени. Для первых опытов оставьте 48 реализаций."),
    h("Как понять состояние задачи"),
    table(["Надпись","Что делать"],[
        ("В очереди","Дождитесь освобождения программы."),
        ("Выполняется / точка ... из ...","Подождите: расчёт продолжается."),
        ("Готово","Нажмите «Открыть результат»."),
        ("Ошибка","Передайте текст сообщения помощнику или коллеге."),
        ("Прервано","Найдите работу в «Результатах» и нажмите «Возобновить», если кнопка доступна.")]),
    Spacer(1,8),p("Пока идёт расчёт, ноутбуком можно пользоваться. Не закрывайте служебное окно запуска и не выключайте устройство. Первый запуск может занять больше времени, чем повторный.")]
story.append(PageBreak())
# 5. Common daily operations; leave advanced data preparation to colleague.
story += [p("Сохранение и свои измерения","title"),h("Сохранить рисунок или таблицу"),
    p("Результаты сохраняются автоматически. Отдельно сохранять расчёт перед выходом не требуется."),
    step(1,"Откройте работу в <b>«Результатах»</b> и раскройте <b>«Таблицы и файлы результата»</b>."),
    step(2,"Для рисунка откройте файл <b>.png</b> и сохраните его обычной командой браузера. Для таблицы используйте файл <b>.csv</b>."),
    p("<b>fmr_vs_field.png</b> - график частоты от поля;<br/><b>fmr_vs_field.csv</b> - таблица к нему;<br/><b>fmr_spectra.png</b> - обзор спектров."),
    p("Если Excel поместил всё в одну колонку: «Данные → Из текста/CSV», разделитель «запятая». Последние работы находятся сверху списка. Флажок «Показывать предыдущие версии» открывает историю; файлы при этом не удаляются."),
    h("Обработать файл измерения"),
    p("<b>1.</b> После прогноза откройте «Проверка модели», выберите его и нажмите «Зафиксировать этот прогноз».<br/><b>2.</b> Перейдите в «Анализ данных» и выберите метод.<br/><b>3.</b> Вставьте полный путь к исходному файлу. В Проводнике: правая кнопка мыши → «Копировать как путь». Кавычки удалять не обязательно.<br/><b>4.</b> Заполните поля выбранного метода и нажмите «Обработать».<br/><b>5.</b> Дождитесь «Готово» и откройте результат."),
    table(["Метод","Что подготовить"],[
        ("РФА","Файл измерения и длину волны."),
        ("Мёссбауэр","Файл спектра; для физических полей и сдвигов - калибровку α-Fe. Без неё показывается график по каналам."),
        ("РЭМ","Изображение, масштаб и область без надписей и линейки."),
        ("EDS","Количественно обработанную таблицу прибора."),
        ("ФМР","Таблицу с осью и сигналом, единицы и тип сигнала.")]),
    Spacer(1,7),p("Первое сравнение измерения с прогнозом подготовьте вместе с коллегой: нужен отдельный файл с измеренными значениями и условиями. Обработка кривой сама по себе ещё не является проверкой точности модели.","small")]
story.append(PageBreak())
# 6. Help and physical folder map.
story += [p("Помощь и расположение файлов","title"),
    table(["Ситуация","Что сделать"],[
        ("Страница не открывается","Запустите START_HEXW.cmd и откройте http://127.0.0.1:8765/."),
        ("Порт занят / WinError 10048","Если HexW открывается, программа уже работает. Закройте окно повторного запуска и пользуйтесь страницей."),
        ("На странице старый экран","Обновите вкладку после исправлений программы."),
        ("Расчёт не начинается","Проверьте поля. Для нескольких точек конец диапазона должен быть больше начала."),
        ("Слишком долго считается","В следующий раз сократите число точек. Для знакомства оставьте готовый пример."),
        ("Ноутбук выключился","Запустите программу снова. Проверьте «Результаты»; для прерванной работы может быть доступно «Возобновить»."),
        ("Непонятная ошибка","Скопируйте текст или сделайте снимок экрана для помощника/коллеги.")]),
    Spacer(1,8),h("Где что лежит"),p("Общая папка: <b>C:\\MIPT\\working_stuff\\Hexx</b>"),
    table(["Что найти","Где искать"],[
        ("Запуск программы","START_HEXW.cmd"),
        ("Эта PDF-инструкция","GUIDE → Инструкция_HexW.pdf"),
        ("Графики и таблицы","runs; удобнее открывать через «Результаты» на сайте."),
        ("Исходные измерения","vault → lab_validation"),
        ("Статьи и исходные материалы","literature; raw → theory_side"),
        ("Научные возможности и ограничения","GUIDE → MODEL_STATUS_V02_RU.md")]),
    Spacer(1,7),p("Папки <b>project, registry, tools, audit</b> - служебные. Для обычной работы изменять их не требуется."),
    p("Чтобы сделать резервную копию, дождитесь завершения задач и скопируйте всю папку Hexx на выбранный носитель. При переносе на другой компьютер для первого запуска может потребоваться помощь с настройкой."),
    box("<b>Завершение работы:</b> дождитесь «Готово», затем можно закрыть вкладку. При обычном самостоятельном запуске можно также закрыть открытое вами служебное окно. Сохранённые результаты останутся на месте.")]

doc=SimpleDocTemplate(str(PDF),pagesize=A4,rightMargin=42,leftMargin=42,topMargin=38,bottomMargin=45,
    title="HexW: инструкция пользователя",author="HexW project",subject="Запуск, первый расчёт Co-Zn, результаты и помощь")
doc.build(story,canvasmaker=PageNumbers)
reader=PdfReader(str(PDF))
texts=[page.extract_text() or "" for page in reader.pages]
assert len(reader.pages)==6, f"Expected six pages, got {len(reader.pages)}"
assert all(len(t)>300 for t in texts)
assert "Быстрый" not in texts[0] or "START_HEXW.cmd" in texts[0]
assert "10048" in texts[0] and "10048" in texts[5]
assert all("�" not in t for t in texts)
report={"file":str(PDF),"pages":len(reader.pages),"bytes":PDF.stat().st_size,"page_text_lengths":[len(t) for t in texts]}
(OUT/"checks.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
print(json.dumps(report,ensure_ascii=False),flush=True)
