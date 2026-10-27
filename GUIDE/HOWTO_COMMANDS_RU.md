# Справочник команд

Программа вызывается как `python -m hexw <команда>`. Все примеры — для PowerShell на Windows.

## 1. Подготовка

Окружение (Python 3.13 или новее) уже создано: `D:\Work\Hexx\.venv`. Пакет `hexw` установлен в режиме «редактируемый» из `D:\Work\Hexx\project`.

    D:\Work\Hexx\.venv\Scripts\Activate.ps1
    cd D:\Work\Hexx\project

Если окружение нужно создать заново (например на другом компьютере; перед этим уточните у разработчика, это не входит в обычную работу):

    python -m venv D:\Work\Hexx\.venv
    D:\Work\Hexx\.venv\Scripts\Activate.ps1
    pip install -r D:\Work\Hexx\project\requirements.txt

В `requirements.txt` пакет `hexw` указан как редактируемая установка по пути `d:\work\hexx\project`; при другом расположении путь нужно поправить.

## 2. Работающие команды

| Команда | Что делает | Код возврата |
|---|---|---|
| `python -m hexw audit` | Проверка целостности и защиты. Выводит `[OK]` или `[FAIL]` по каждой проверке и итог. Ничего не записывает. | 0 — всё хорошо; 1 — есть провалы |
| `python -m hexw audit --write-manifest` | То же и сохраняет манифест запуска (хэши, версии git, окружение) в `D:\Work\Hexx\registry\runs\` | 0 / 1 |
| `python -m hexw audit --root <путь>` | Явно указать корень проекта; он обязан совпадать с настроенным (иначе задаётся переменной `HEXW_ROOT`). Обычно не нужен. | 0 / 1 / 2 |
| `python -m pytest` | Автотесты проекта (из папки `project`). Нагружает систему — не запускать одновременно с тяжёлыми расчётами. | стандарт pytest |

Справка: `python -m hexw --help` и `python -m hexw audit --help`.

Аудит и тесты полезно запускать после любых изменений в проекте и перед важными этапами (фиксация модели, сравнение с данными).

## 3. Команды-заглушки

Они присутствуют, но возвращают «этап не реализован. Результатов нет.» и код 3:

    python -m hexw literature      # литературный корпус (этап 2)
    python -m hexw structures      # структурные кандидаты (этап 4)
    python -m hexw dft             # DFT-параметризация (этап 5)
    python -m hexw spin-model      # спиновая модель / ФМР (этапы 4–6)
    python -m hexw predict         # предварительные прогнозы для серии (этап 6)
    python -m hexw freeze-model    # фиксация модели и прогнозов (этап 7)
    python -m hexw analyze-exp     # независимый анализ экспериментов (этап 8)
    python -m hexw validate        # сравнение теории и эксперимента (этап 9)
    python -m hexw report          # итоговый отчёт
    python -m hexw resume          # возобновление с контрольной точки

Номера этапов в справке программы соответствуют разделу 15 ТЗ.

**Калибровка по лабораторным данным:**

    python -m hexw explicit-calibrate

без подтверждения отказывает (код 2) и ничего не делает. Подтверждающий флаг: `--i-confirm-explicit-calibration`. Сейчас и с флагом это заглушка. Используйте только осознанно:
калибровка создаёт новую версию модели и лишает использованные данные статуса независимой проверки.

## 4. Что где смотреть

| Что | Где |
|---|---|
| Отчёты аудита | `D:\Work\Hexx\registry\runs\` |
| Образцы, измерения, прежнее использование | `D:\Work\Hexx\registry\samples.yaml`, `observations.yaml`, `prior_use.yaml` |
| Литературные записи для подгонки | `D:\Work\Hexx\literature\fit\records_fit.yaml` |
| Пробелы литературы | `D:\Work\Hexx\literature\fit\coverage_map.csv` |
| Противоречия источников | `D:\Work\Hexx\literature\fit\CONFLICTS.md` |
| Нотация позиций W1...W7 | `D:\Work\Hexx\literature\fit\site_notation.md` |
| Решение о разбиении fit / holdout | `D:\Work\Hexx\literature\SPLIT_DECISION.md` |
| Журнал литературного поиска и кандидаты | `D:\Work\Hexx\lit_search\SEARCH_LOG.md`, `lit_search\results\candidates.csv` |
| Возможности DFT на ноутбуке | `D:\Work\Hexx\dft_feasibility\FEASIBILITY.md`, `install_log.md` |
| Описание изоляции данных | `D:\Work\Hexx\project\docs\ISOLATION.md` |
| Описание магнитной модели | `D:\Work\Hexx\project\docs\MAGNET_MODEL.md` |

## 5. Что нельзя делать

- Не править и не удалять содержимое `raw\`, `vault\`, `literature_holdout\` и ключевые файлы реестра.
- Не копировать лабораторные данные в папки теории (`project\hexw\theory` и др.): код теории всё равно получит отказ.
- Не запускать скрипты разбиения литературы с ключом `--force` без согласования: исходные смешанные записи лежат в защищённом резерве.
- Не запускать `explicit-calibrate` «на пробу».

## 6. Если аудит сообщает о провале

1. Прочитайте строку `[FAIL] <проверка> -> <детали>`.
2. Не пытайтесь «починить» правкой защищённых файлов. Если файл изменён случайно — верните его из git (`registry`) или из исходного архива в `raw\original_zip`.
3. Обратитесь к разработчику, приложив вывод команды.
