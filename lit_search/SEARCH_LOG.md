# SEARCH_LOG

Дата выполнения: 2026-10-25 (UTC). User-Agent: `hexw-litsearch/0.1`, без ключей, без e-mail, без личных данных в запросах. Unpaywall не использовался.

## Сетевая проверка
`curl https://api.openalex.org/works?search=W-type hexaferrite&per-page=1` -> HTTP 200, JSON получен. Интернет есть.

## Ограничения и инциденты
- OpenAlex: без ключа общий дневной бюджет ($0.10 на IP) был исчерпан после первых ~22 запросов (по теме 1, 2 и части 3); остальные запросы вернули HTTP 429 (в таблице ниже -- 0 результатов, это отказ, а не пустая выдача). Ключ не запрашивали (требует e-mail пользователя). Эти темы покрыты Crossref, Europe PMC, DOAJ, arXiv и частично Semantic Scholar.
- Semantic Scholar graph API: 429 почти на всех запросах (из 29 поисковых -- 1 успешный; паузы 6-8 с и 2 повтора не помогли). Использовать как источник нельзя; для L04/L06/L07/L08 точечные запросы по DOI прошли (openAccessPdf=CLOSED).
- arXiv API: фразовые запросы в кавычках для длинных строк часто дают 0; добавлен проход простыми запросами (all:hexaferrite -> 112 записей и др.).
- Crossref `query.bibliographic` и `query` -- шумные (total-results в миллионах), использовались только как источник кандидатов с последующей локальной фильтрацией по заголовку/реферату.
- Рефераты есть не у всех записей (Elsevier/Springer/AIP часто без реферата в Crossref); `abstract_available=no` -- оценка релевантности только по заголовку.
- Сайты MDPI, T&F, RSC и europepmc render-endpoint отвечают 403 на нейтральный бот-UA; обход (смена UA, прокси) не применялся. IOP -- SSL EOF. Часть OA-статей поэтому не скачана; получены через Europe PMC REST fullTextXML, где статья в PMC OA.
- J-STAGE: 3 PDF (C064, C069, C070) получены по OA-ссылке, но не верифицируются автоматически (скан/японские шрифты, текст не извлекается), поэтому удалены.
- OA-статус определялся по метаданным источников (OpenAlex is_oa/oa_url, Europe PMC isOpenAccess, DOAJ, arXiv, Semantic Scholar openAccessPdf) и по CC-лицензии в Crossref. Колонка oa_status отражает метаданные, а не проверку условий лицензии.
- Ни один платный текст не получен; Sci-Hub и подобные не использовались.

## Запросы (источник, запрос, число результатов по данным API на момент запроса)

| файл | API | запрос | n | время UTC |
|---|---|---|---|---|
| pass | OpenAlex | SrCo2W hexaferrite | 4 | 14:41:04 |
| pass | Crossref | SrCo2W hexaferrite | 2850 | 14:41:04 |
| pass | OpenAlex | SrZn2W hexaferrite | 3 | 14:41:04 |
| pass | Crossref | SrZn2W hexaferrite | 2852 | 14:41:04 |
| pass | OpenAlex | SrMn2W hexaferrite | 0 | 14:41:04 |
| pass | Crossref | SrMn2W hexaferrite | 2850 | 14:41:04 |
| pass | OpenAlex | SrFe18O27 | 19 | 14:41:04 |
| pass | Crossref | SrFe18O27 | 2 | 14:41:04 |
| pass | OpenAlex | Sr W-type hexaferrite | 2063 | 14:41:04 |
| pass | Crossref | Sr W-type hexaferrite | 1383209 | 14:41:04 |
| pass | OpenAlex | Sr(Co,Zn)2Fe16O27 | 9 | 14:41:04 |
| pass | Crossref | Sr(Co,Zn)2Fe16O27 | 1094473 | 14:41:04 |
| pass | OpenAlex | SrMe2Fe16O27 W-type | 3 | 14:41:04 |
| pass | Crossref | SrMe2Fe16O27 W-type | 1300277 | 14:41:04 |
| pass | OpenAlex | strontium W-type ferrite cobalt zinc substitution | 1141 | 14:41:04 |
| pass | Crossref | strontium W-type ferrite cobalt zinc substitution | 1782886 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite Ni2W | 2 | 14:41:04 |
| pass | Crossref | W-type hexaferrite Ni2W | 1302428 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite Mg2W | 2 | 14:41:04 |
| pass | Crossref | W-type hexaferrite Mg2W | 1302428 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite Cu2W | 1 | 14:41:04 |
| pass | Crossref | W-type hexaferrite Cu2W | 1302428 | 14:41:04 |
| pass | OpenAlex | Fe2W hexaferrite | 8 | 14:41:04 |
| pass | Crossref | Fe2W hexaferrite | 2870 | 14:41:04 |
| pass | OpenAlex | W-type ferrite BaFe2+2Fe16O27 substitution | 2 | 14:41:04 |
| pass | Crossref | W-type ferrite BaFe2+2Fe16O27 substitution | 1427696 | 14:41:04 |
| pass | OpenAlex | BaCo2W hexaferrite | 3 | 14:41:04 |
| pass | Crossref | BaCo2W hexaferrite | 2850 | 14:41:04 |
| pass | OpenAlex | BaZn2W hexaferrite | 2 | 14:41:04 |
| pass | Crossref | BaZn2W hexaferrite | 2850 | 14:41:04 |
| pass | OpenAlex | BaMn2W hexaferrite | 0 | 14:41:04 |
| pass | Crossref | BaMn2W hexaferrite | 2850 | 14:41:04 |
| pass | OpenAlex | BaFe18O27 | 29 | 14:41:04 |
| pass | Crossref | BaFe18O27 | 9 | 14:41:04 |
| pass | OpenAlex | Co-Zn W-type hexaferrite | 2048 | 14:41:04 |
| pass | Crossref | Co-Zn W-type hexaferrite | 2306099 | 14:41:04 |
| pass | OpenAlex | spin reorientation W-type hexaferrite | 194 | 14:41:04 |
| pass | Crossref | spin reorientation W-type hexaferrite | 1561674 | 14:41:04 |
| pass | OpenAlex | cone anisotropy W-type hexaferrite | 139 | 14:41:04 |
| pass | Crossref | cone anisotropy W-type hexaferrite | 1435310 | 14:41:04 |
| pass | OpenAlex | ferromagnetic resonance W-type hexaferrite | 1034 | 14:41:04 |
| pass | Crossref | ferromagnetic resonance W-type hexaferrite | 1789351 | 14:41:04 |
| pass | OpenAlex | microwave permeability W-type hexaferrite natural resonance | 432 | 14:41:04 |
| pass | Crossref | microwave permeability W-type hexaferrite natural resonance | 3249130 | 14:41:04 |
| pass | OpenAlex | Mossbauer W-type hexaferrite site occupancy | 0 | 14:41:04 |
| pass | Crossref | Mossbauer W-type hexaferrite site occupancy | 1687606 | 14:41:04 |
| pass | OpenAlex | cation distribution W-type hexaferrite | 0 | 14:41:04 |
| pass | Crossref | cation distribution W-type hexaferrite | 2257396 | 14:41:04 |
| pass | OpenAlex | neutron diffraction W-type hexaferrite | 0 | 14:41:04 |
| pass | Crossref | neutron diffraction W-type hexaferrite | 1572244 | 14:41:04 |
| pass | OpenAlex | first-principles SrFe18O27 | 0 | 14:41:04 |
| pass | Crossref | first-principles SrFe18O27 | 1485023 | 14:41:04 |
| pass | OpenAlex | first-principles BaFe18O27 W-type | 0 | 14:41:04 |
| pass | Crossref | first-principles BaFe18O27 W-type | 2774204 | 14:41:04 |
| pass | OpenAlex | exchange constants hexaferrite DFT+U | 0 | 14:41:04 |
| pass | Crossref | exchange constants hexaferrite DFT+U | 527297 | 14:41:04 |
| pass | OpenAlex | magnetocrystalline anisotropy hexaferrite spin-orbit first principles | 0 | 14:41:04 |
| pass | Crossref | magnetocrystalline anisotropy hexaferrite spin-orbit first principles | 1850997 | 14:41:04 |
| pass | OpenAlex | SrFe12O19 density functional theory magnetic anisotropy | 0 | 14:41:04 |
| pass | Crossref | SrFe12O19 density functional theory magnetic anisotropy | 4268598 | 14:41:04 |
| pass | OpenAlex | M-type hexaferrite exchange interactions first-principles | 0 | 14:41:04 |
| pass | Crossref | M-type hexaferrite exchange interactions first-principles | 3700279 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite dielectric permittivity | 0 | 14:41:04 |
| pass | Crossref | W-type hexaferrite dielectric permittivity | 1507952 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 0 | 14:41:04 |
| pass | Crossref | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 2379229 | 14:41:04 |
| pass | OpenAlex | W-type hexaferrite dielectric properties | 0 | 14:41:04 |
| pass | Crossref | W-type hexaferrite dielectric properties | 3498184 | 14:41:04 |
| pass | OpenAlex | Mossbauer recoil-free fraction site populations ferrites | 0 | 14:41:04 |
| pass | Crossref | Mossbauer recoil-free fraction site populations ferrites | 1520320 | 14:41:04 |
| pass | OpenAlex | Smit-Beljers ferromagnetic resonance hexagonal anisotropy | 0 | 14:41:04 |
| pass | Crossref | Smit-Beljers ferromagnetic resonance hexagonal anisotropy | 619847 | 14:41:04 |
| pass | OpenAlex | polycrystalline ferromagnetic resonance linewidth Schlomann porosity | 0 | 14:41:04 |
| pass | Crossref | polycrystalline ferromagnetic resonance linewidth Schlomann porosity | 569767 | 14:41:04 |
| pass | OpenAlex | ferromagnetic resonance hexagonal ferrite single crystal anisotropy field | 0 | 14:41:04 |
| pass | Crossref | ferromagnetic resonance hexagonal ferrite single crystal anisotropy field | 3522860 | 14:41:04 |
| pass2 | EuropePMC | SrCo2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | SrCo2W hexaferrite | 1 | 14:45:35 |
| pass2 | Crossref-query | SrCo2W hexaferrite | 2502 | 14:45:35 |
| pass2 | DOAJ | SrCo2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | SrZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | SrZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | SrZn2W hexaferrite | 2504 | 14:45:35 |
| pass2 | DOAJ | SrZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | SrMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | SrMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | SrMn2W hexaferrite | 2502 | 14:45:35 |
| pass2 | DOAJ | SrMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | SrFe18O27 | 1 | 14:45:35 |
| pass2 | arXiv | SrFe18O27 | 0 | 14:45:35 |
| pass2 | Crossref-query | SrFe18O27 | 2 | 14:45:35 |
| pass2 | DOAJ | SrFe18O27 | 1 | 14:45:35 |
| pass2 | EuropePMC | Sr W-type hexaferrite | 74 | 14:45:35 |
| pass2 | arXiv | Sr W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | Sr W-type hexaferrite | 1182252 | 14:45:35 |
| pass2 | DOAJ | Sr W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | Sr(Co,Zn)2Fe16O27 | 0 | 14:45:35 |
| pass2 | arXiv | Sr(Co,Zn)2Fe16O27 | 12290 | 14:45:35 |
| pass2 | Crossref-query | Sr(Co,Zn)2Fe16O27 | 1065396 | 14:45:35 |
| pass2 | DOAJ | Sr(Co,Zn)2Fe16O27 | 0 | 14:45:35 |
| pass2 | EuropePMC | SrMe2Fe16O27 W-type | 1 | 14:45:35 |
| pass2 | arXiv | SrMe2Fe16O27 W-type | 0 | 14:45:35 |
| pass2 | Crossref-query | SrMe2Fe16O27 W-type | 1071932 | 14:45:35 |
| pass2 | DOAJ | SrMe2Fe16O27 W-type | 0 | 14:45:35 |
| pass2 | EuropePMC | strontium W-type ferrite cobalt zinc substitution | 21 | 14:45:35 |
| pass2 | arXiv | strontium W-type ferrite cobalt zinc substitution | 0 | 14:45:35 |
| pass2 | Crossref-query | strontium W-type ferrite cobalt zinc substitution | 1478341 | 14:45:35 |
| pass2 | DOAJ | strontium W-type ferrite cobalt zinc substitution | 0 | 14:45:35 |
| pass2 | EuropePMC | Sr W-type hexaferrite microwave absorption | 34 | 14:45:35 |
| pass2 | arXiv | Sr W-type hexaferrite microwave absorption | 0 | 14:45:35 |
| pass2 | Crossref-query | Sr W-type hexaferrite microwave absorption | 1633976 | 14:45:35 |
| pass2 | DOAJ | Sr W-type hexaferrite microwave absorption | 0 | 14:45:35 |
| pass2 | EuropePMC | SrCo2W ferrite | 0 | 14:45:35 |
| pass2 | arXiv | SrCo2W ferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | SrCo2W ferrite | 29903 | 14:45:35 |
| pass2 | DOAJ | SrCo2W ferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | Zn2W ferrite | 0 | 14:45:35 |
| pass2 | arXiv | Zn2W ferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | Zn2W ferrite | 29910 | 14:45:35 |
| pass2 | DOAJ | Zn2W ferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite Ni2W | 0 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite Ni2W | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite Ni2W | 1073801 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite Ni2W | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite Mg2W | 0 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite Mg2W | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite Mg2W | 1073801 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite Mg2W | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite Cu2W | 0 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite Cu2W | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite Cu2W | 1073801 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite Cu2W | 0 | 14:45:35 |
| pass2 | EuropePMC | Fe2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | Fe2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | Fe2W hexaferrite | 2522 | 14:45:35 |
| pass2 | DOAJ | Fe2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type ferrite BaFe2+2Fe16O27 substitution | 0 | 14:45:35 |
| pass2 | arXiv | W-type ferrite BaFe2+2Fe16O27 substitution | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type ferrite BaFe2+2Fe16O27 substitution | 1176350 | 14:45:35 |
| pass2 | DOAJ | W-type ferrite BaFe2+2Fe16O27 substitution | 0 | 14:45:35 |
| pass2 | EuropePMC | BaCo2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | BaCo2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | BaCo2W hexaferrite | 2502 | 14:45:35 |
| pass2 | DOAJ | BaCo2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | BaZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | BaZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | BaZn2W hexaferrite | 2502 | 14:45:35 |
| pass2 | DOAJ | BaZn2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | BaMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | arXiv | BaMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | BaMn2W hexaferrite | 2502 | 14:45:35 |
| pass2 | DOAJ | BaMn2W hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | BaFe18O27 | 2 | 14:45:35 |
| pass2 | arXiv | BaFe18O27 | 0 | 14:45:35 |
| pass2 | Crossref-query | BaFe18O27 | 8 | 14:45:35 |
| pass2 | DOAJ | BaFe18O27 | 0 | 14:45:35 |
| pass2 | EuropePMC | Co-Zn W-type hexaferrite | 43 | 14:45:35 |
| pass2 | arXiv | Co-Zn W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | Co-Zn W-type hexaferrite | 2022477 | 14:45:35 |
| pass2 | DOAJ | Co-Zn W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | spin reorientation W-type hexaferrite | 11 | 14:45:35 |
| pass2 | arXiv | spin reorientation W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | spin reorientation W-type hexaferrite | 1290117 | 14:45:35 |
| pass2 | DOAJ | spin reorientation W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | cone anisotropy W-type hexaferrite | 7 | 14:45:35 |
| pass2 | arXiv | cone anisotropy W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | cone anisotropy W-type hexaferrite | 1177848 | 14:45:35 |
| pass2 | DOAJ | cone anisotropy W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | ferromagnetic resonance W-type hexaferrite | 40 | 14:45:35 |
| pass2 | arXiv | ferromagnetic resonance W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | ferromagnetic resonance W-type hexaferrite | 1479453 | 14:45:35 |
| pass2 | DOAJ | ferromagnetic resonance W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | microwave permeability W-type hexaferrite natural resonance | 22 | 14:45:35 |
| pass2 | arXiv | microwave permeability W-type hexaferrite natural resonance | 0 | 14:45:35 |
| pass2 | Crossref-query | microwave permeability W-type hexaferrite natural resonance | 6349439 | 14:45:35 |
| pass2 | DOAJ | microwave permeability W-type hexaferrite natural resonance | 0 | 14:45:35 |
| pass2 | EuropePMC | Mossbauer W-type hexaferrite site occupancy | 6 | 14:45:35 |
| pass2 | arXiv | Mossbauer W-type hexaferrite site occupancy | 0 | 14:45:35 |
| pass2 | Crossref-query | Mossbauer W-type hexaferrite site occupancy | 1352312 | 14:45:35 |
| pass2 | DOAJ | Mossbauer W-type hexaferrite site occupancy | 0 | 14:45:35 |
| pass2 | EuropePMC | cation distribution W-type hexaferrite | 57 | 14:45:35 |
| pass2 | arXiv | cation distribution W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | cation distribution W-type hexaferrite | 1791878 | 14:45:35 |
| pass2 | DOAJ | cation distribution W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | neutron diffraction W-type hexaferrite | 16 | 14:45:35 |
| pass2 | arXiv | neutron diffraction W-type hexaferrite | 0 | 14:45:35 |
| pass2 | Crossref-query | neutron diffraction W-type hexaferrite | 1286875 | 14:45:35 |
| pass2 | DOAJ | neutron diffraction W-type hexaferrite | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite Mossbauer | 12 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite Mossbauer | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite Mossbauer | 1075024 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite Mossbauer | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite magnetic structure | 124 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite magnetic structure | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite magnetic structure | 3281311 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite magnetic structure | 0 | 14:45:35 |
| pass2 | EuropePMC | barium W-type ferrite single crystal magnetic anisotropy | 92 | 14:45:35 |
| pass2 | arXiv | barium W-type ferrite single crystal magnetic anisotropy | 0 | 14:45:35 |
| pass2 | Crossref-query | barium W-type ferrite single crystal magnetic anisotropy | 3214719 | 14:45:35 |
| pass2 | DOAJ | barium W-type ferrite single crystal magnetic anisotropy | 0 | 14:45:35 |
| pass2 | EuropePMC | hexagonal ferrite Co2W anisotropy temperature dependence | 0 | 14:45:35 |
| pass2 | arXiv | hexagonal ferrite Co2W anisotropy temperature dependence | 0 | 14:45:35 |
| pass2 | Crossref-query | hexagonal ferrite Co2W anisotropy temperature dependence | 1106792 | 14:45:35 |
| pass2 | DOAJ | hexagonal ferrite Co2W anisotropy temperature dependence | 0 | 14:45:35 |
| pass2 | EuropePMC | first-principles SrFe18O27 | 0 | 14:45:35 |
| pass2 | arXiv | first-principles SrFe18O27 | 0 | 14:45:35 |
| pass2 | Crossref-query | first-principles SrFe18O27 | 889000 | 14:45:35 |
| pass2 | EuropePMC | first-principles BaFe18O27 W-type | 0 | 14:45:35 |
| pass2 | arXiv | first-principles BaFe18O27 W-type | 0 | 14:45:35 |
| pass2 | Crossref-query | first-principles BaFe18O27 W-type | 1952437 | 14:45:35 |
| pass2 | EuropePMC | exchange constants hexaferrite DFT+U | 9 | 14:45:35 |
| pass2 | arXiv | exchange constants hexaferrite DFT+U | 0 | 14:45:35 |
| pass2 | Crossref-query | exchange constants hexaferrite DFT+U | 395071 | 14:45:35 |
| pass2 | EuropePMC | magnetocrystalline anisotropy hexaferrite spin-orbit first principles | 5 | 14:45:35 |
| pass2 | arXiv | magnetocrystalline anisotropy hexaferrite spin-orbit first principles | 0 | 14:45:35 |
| pass2 | Crossref-query | magnetocrystalline anisotropy hexaferrite spin-orbit first principles | 1183624 | 14:45:35 |
| pass2 | EuropePMC | SrFe12O19 density functional theory magnetic anisotropy | 20 | 14:45:35 |
| pass2 | arXiv | SrFe12O19 density functional theory magnetic anisotropy | 0 | 14:45:35 |
| pass2 | Crossref-query | SrFe12O19 density functional theory magnetic anisotropy | 2997912 | 14:45:35 |
| pass2 | EuropePMC | M-type hexaferrite exchange interactions first-principles | 24 | 14:45:35 |
| pass2 | arXiv | M-type hexaferrite exchange interactions first-principles | 0 | 14:45:35 |
| pass2 | Crossref-query | M-type hexaferrite exchange interactions first-principles | 2653700 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite dielectric permittivity | 36 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite dielectric permittivity | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite dielectric permittivity | 1220040 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite dielectric permittivity | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 13 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 1693157 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite electrical conductivity Fe2+ Fe3+ hopping | 0 | 14:45:35 |
| pass2 | EuropePMC | W-type hexaferrite dielectric properties | 64 | 14:45:35 |
| pass2 | arXiv | W-type hexaferrite dielectric properties | 0 | 14:45:35 |
| pass2 | Crossref-query | W-type hexaferrite dielectric properties | 2796524 | 14:45:35 |
| pass2 | DOAJ | W-type hexaferrite dielectric properties | 0 | 14:45:35 |
| pass2 | EuropePMC | Mossbauer recoil-free fraction site populations ferrites | 3 | 14:45:35 |
| pass2 | arXiv | Mossbauer recoil-free fraction site populations ferrites | 0 | 14:45:35 |
| pass2 | Crossref-query | Mossbauer recoil-free fraction site populations ferrites | 1158140 | 14:45:35 |
| pass2 | EuropePMC | Smit-Beljers ferromagnetic resonance hexagonal anisotropy | 0 | 14:45:35 |
| pass2 | arXiv | Smit-Beljers ferromagnetic resonance hexagonal anisotropy | 0 | 14:45:35 |
| pass2 | Crossref-query | Smit-Beljers ferromagnetic resonance hexagonal anisotropy | 510287 | 14:45:35 |
| pass2 | EuropePMC | polycrystalline ferromagnetic resonance linewidth Schlomann porosity | 0 | 14:45:35 |
| pass2 | arXiv | polycrystalline ferromagnetic resonance linewidth Schlomann porosity | 0 | 14:45:35 |
| pass2 | Crossref-query | polycrystalline ferromagnetic resonance linewidth Schlomann porosity | 471946 | 14:45:35 |
| pass2 | EuropePMC | ferromagnetic resonance hexagonal ferrite single crystal anisotropy field | 125 | 14:45:35 |
| pass2 | arXiv | ferromagnetic resonance hexagonal ferrite single crystal anisotropy field | 0 | 14:45:35 |
| pass2 | Crossref-query | ferromagnetic resonance hexagonal ferrite single crystal anisotropy field | 2643821 | 14:45:35 |
| pass3 | EuropePMC | ab initio magnetic anisotropy strontium hexaferrite | 3 | 14:51:30 |
| pass3 | arXiv | ab initio magnetic anisotropy strontium hexaferrite | 0 | 14:51:30 |
| pass3 | Crossref-query | ab initio magnetic anisotropy strontium hexaferrite | 962875 | 14:51:30 |
| pass3 | SemanticScholar | ab initio magnetic anisotropy strontium hexaferrite | 0 | 14:51:30 |
| pass3 | EuropePMC | SrFe12O19 first principles electronic structure | 32 | 14:51:30 |
| pass3 | arXiv | SrFe12O19 first principles electronic structure | 0 | 14:51:30 |
| pass3 | Crossref-query | SrFe12O19 first principles electronic structure | 3054172 | 14:51:30 |
| pass3 | SemanticScholar | SrFe12O19 first principles electronic structure | 0 | 14:51:30 |
| pass3 | EuropePMC | BaFe12O19 DFT+U exchange interactions | 9 | 14:51:30 |
| pass3 | arXiv | BaFe12O19 DFT+U exchange interactions | 0 | 14:51:30 |
| pass3 | Crossref-query | BaFe12O19 DFT+U exchange interactions | 762610 | 14:51:30 |
| pass3 | SemanticScholar | BaFe12O19 DFT+U exchange interactions | 0 | 14:51:30 |
| pass3 | EuropePMC | hexaferrite magnetocrystalline anisotropy density functional theory | 20 | 14:51:30 |
| pass3 | arXiv | hexaferrite magnetocrystalline anisotropy density functional theory | 0 | 14:51:30 |
| pass3 | Crossref-query | hexaferrite magnetocrystalline anisotropy density functional theory | 2235562 | 14:51:30 |
| pass3 | SemanticScholar | hexaferrite magnetocrystalline anisotropy density functional theory | 5747 | 14:51:30 |
| pass3 | EuropePMC | M-type hexaferrite substitution first-principles Co Zn | 7 | 14:51:30 |
| pass3 | arXiv | M-type hexaferrite substitution first-principles Co Zn | 0 | 14:51:30 |
| pass3 | Crossref-query | M-type hexaferrite substitution first-principles Co Zn | 2965757 | 14:51:30 |
| pass3 | SemanticScholar | M-type hexaferrite substitution first-principles Co Zn | 0 | 14:51:30 |
| pass3 | EuropePMC | W-type hexaferrite first-principles | 22 | 14:51:30 |
| pass3 | arXiv | W-type hexaferrite first-principles | 0 | 14:51:30 |
| pass3 | Crossref-query | W-type hexaferrite first-principles | 1954293 | 14:51:30 |
| pass3 | SemanticScholar | W-type hexaferrite first-principles | 0 | 14:51:30 |
| pass3 | EuropePMC | Heisenberg exchange Monte Carlo hexaferrite Curie temperature | 0 | 14:51:30 |
| pass3 | arXiv | Heisenberg exchange Monte Carlo hexaferrite Curie temperature | 0 | 14:51:30 |
| pass3 | Crossref-query | Heisenberg exchange Monte Carlo hexaferrite Curie temperature | 1456192 | 14:51:30 |
| pass3 | SemanticScholar | Heisenberg exchange Monte Carlo hexaferrite Curie temperature | 0 | 14:51:30 |
| pass3 | EuropePMC | hexaferrite spin-orbit coupling magnetic anisotropy site-resolved | 7 | 14:51:30 |
| pass3 | arXiv | hexaferrite spin-orbit coupling magnetic anisotropy site-resolved | 0 | 14:51:30 |
| pass3 | Crossref-query | hexaferrite spin-orbit coupling magnetic anisotropy site-resolved | 1592897 | 14:51:30 |
| pass3 | SemanticScholar | hexaferrite spin-orbit coupling magnetic anisotropy site-resolved | 0 | 14:51:30 |
| pass3 | EuropePMC | W-type hexaferrite magnetic structure neutron | 18 | 14:51:30 |
| pass3 | arXiv | W-type hexaferrite magnetic structure neutron | 0 | 14:51:30 |
| pass3 | Crossref-query | W-type hexaferrite magnetic structure neutron | 3395130 | 14:51:30 |
| pass3 | DOAJ | W-type hexaferrite magnetic structure neutron | 0 | 14:51:30 |
| pass3 | SemanticScholar | W-type hexaferrite magnetic structure neutron | 0 | 14:51:30 |
| pass3 | EuropePMC | Co2W single crystal anisotropy | 2 | 14:51:30 |
| pass3 | arXiv | Co2W single crystal anisotropy | 0 | 14:51:30 |
| pass3 | Crossref-query | Co2W single crystal anisotropy | 1361557 | 14:51:30 |
| pass3 | DOAJ | Co2W single crystal anisotropy | 0 | 14:51:30 |
| pass3 | SemanticScholar | Co2W single crystal anisotropy | 0 | 14:51:30 |
| pass3 | EuropePMC | Fe2+ site occupancy W-type Mossbauer | 140 | 14:51:30 |
| pass3 | arXiv | Fe2+ site occupancy W-type Mossbauer | 0 | 14:51:30 |
| pass3 | Crossref-query | Fe2+ site occupancy W-type Mossbauer | 1354967 | 14:51:30 |
| pass3 | DOAJ | Fe2+ site occupancy W-type Mossbauer | 0 | 14:51:30 |
| pass3 | SemanticScholar | Fe2+ site occupancy W-type Mossbauer | 0 | 14:51:30 |
| pass3 | EuropePMC | W-type hexaferrite thin film ferromagnetic resonance | 11 | 14:51:30 |
| pass3 | arXiv | W-type hexaferrite thin film ferromagnetic resonance | 0 | 14:51:30 |
| pass3 | Crossref-query | W-type hexaferrite thin film ferromagnetic resonance | 2122044 | 14:51:30 |
| pass3 | DOAJ | W-type hexaferrite thin film ferromagnetic resonance | 0 | 14:51:30 |
| pass3 | SemanticScholar | W-type hexaferrite thin film ferromagnetic resonance | 0 | 14:51:30 |
| pass3 | EuropePMC | self-biased circulator hexaferrite W-type | 1 | 14:51:30 |
| pass3 | arXiv | self-biased circulator hexaferrite W-type | 0 | 14:51:30 |
| pass3 | Crossref-query | self-biased circulator hexaferrite W-type | 1988911 | 14:51:30 |
| pass3 | DOAJ | self-biased circulator hexaferrite W-type | 0 | 14:51:30 |
| pass3 | SemanticScholar | self-biased circulator hexaferrite W-type | 0 | 14:51:30 |
| pass3 | EuropePMC | BaCo2Fe16O27 | 4 | 14:51:30 |
| pass3 | arXiv | BaCo2Fe16O27 | 0 | 14:51:30 |
| pass3 | Crossref-query | BaCo2Fe16O27 | 8 | 14:51:30 |
| pass3 | DOAJ | BaCo2Fe16O27 | 1 | 14:51:30 |
| pass3 | SemanticScholar | BaCo2Fe16O27 | 0 | 14:51:30 |
| pass3 | EuropePMC | BaZnCoFe16O27 magnetic anisotropy | 0 | 14:51:30 |
| pass3 | arXiv | BaZnCoFe16O27 magnetic anisotropy | 0 | 14:51:30 |
| pass3 | Crossref-query | BaZnCoFe16O27 magnetic anisotropy | 843746 | 14:51:30 |
| pass3 | DOAJ | BaZnCoFe16O27 magnetic anisotropy | 0 | 14:51:30 |
| pass3 | SemanticScholar | BaZnCoFe16O27 magnetic anisotropy | 0 | 14:51:30 |
| pass3 | EuropePMC | strontium W-type hexaferrite Mossbauer | 7 | 14:51:30 |
| pass3 | arXiv | strontium W-type hexaferrite Mossbauer | 0 | 14:51:30 |
| pass3 | Crossref-query | strontium W-type hexaferrite Mossbauer | 1097757 | 14:51:30 |
| pass3 | DOAJ | strontium W-type hexaferrite Mossbauer | 0 | 14:51:30 |
| pass3 | SemanticScholar | strontium W-type hexaferrite Mossbauer | 0 | 14:51:30 |
| pass3 | EuropePMC | strontium W-type hexaferrite neutron diffraction | 5 | 14:51:30 |
| pass3 | arXiv | strontium W-type hexaferrite neutron diffraction | 0 | 14:51:30 |
| pass3 | Crossref-query | strontium W-type hexaferrite neutron diffraction | 1309441 | 14:51:30 |
| pass3 | DOAJ | strontium W-type hexaferrite neutron diffraction | 0 | 14:51:30 |
| pass3 | SemanticScholar | strontium W-type hexaferrite neutron diffraction | 0 | 14:51:30 |
| pass3 | EuropePMC | SrCo2Fe16O27 magnetic | 2 | 14:51:30 |
| pass3 | arXiv | SrCo2Fe16O27 magnetic | 0 | 14:51:30 |
| pass3 | Crossref-query | SrCo2Fe16O27 magnetic | 802560 | 14:51:30 |
| pass3 | DOAJ | SrCo2Fe16O27 magnetic | 0 | 14:51:30 |
| pass3 | SemanticScholar | SrCo2Fe16O27 magnetic | 0 | 14:51:30 |
| pass3 | EuropePMC | Zn2W hexaferrite magnetization cobalt | 0 | 14:51:30 |
| pass3 | arXiv | Zn2W hexaferrite magnetization cobalt | 0 | 14:51:30 |
| pass3 | Crossref-query | Zn2W hexaferrite magnetization cobalt | 133776 | 14:51:30 |
| pass3 | DOAJ | Zn2W hexaferrite magnetization cobalt | 0 | 14:51:30 |
| pass3 | SemanticScholar | Zn2W hexaferrite magnetization cobalt | 0 | 14:51:30 |
| pass3 | EuropePMC | W-type ferrite spin reorientation temperature Co content | 41 | 14:51:30 |
| pass3 | arXiv | W-type ferrite spin reorientation temperature Co content | 0 | 14:51:30 |
| pass3 | Crossref-query | W-type ferrite spin reorientation temperature Co content | 3314992 | 14:51:30 |
| pass3 | DOAJ | W-type ferrite spin reorientation temperature Co content | 0 | 14:51:30 |
| pass3 | SemanticScholar | W-type ferrite spin reorientation temperature Co content | 0 | 14:51:30 |
| pass3 | EuropePMC | dielectric relaxation Fe2+ Fe3+ hopping hexaferrite | 19 | 14:51:30 |
| pass3 | arXiv | dielectric relaxation Fe2+ Fe3+ hopping hexaferrite | 0 | 14:51:30 |
| pass3 | Crossref-query | dielectric relaxation Fe2+ Fe3+ hopping hexaferrite | 265635 | 14:51:30 |
| pass3 | DOAJ | dielectric relaxation Fe2+ Fe3+ hopping hexaferrite | 0 | 14:51:30 |
| pass3 | SemanticScholar | dielectric relaxation Fe2+ Fe3+ hopping hexaferrite | 0 | 14:51:30 |
| pass3 | EuropePMC | dielectric W-type barium ferrite | 184 | 14:51:30 |
| pass3 | arXiv | dielectric W-type barium ferrite | 0 | 14:51:30 |
| pass3 | Crossref-query | dielectric W-type barium ferrite | 1264743 | 14:51:30 |
| pass3 | DOAJ | dielectric W-type barium ferrite | 0 | 14:51:30 |
| pass3 | SemanticScholar | dielectric W-type barium ferrite | 0 | 14:51:30 |
| pass3 | EuropePMC | W-type strontium hexaferrite dielectric | 32 | 14:51:30 |
| pass3 | arXiv | W-type strontium hexaferrite dielectric | 0 | 14:51:30 |
| pass3 | Crossref-query | W-type strontium hexaferrite dielectric | 1236014 | 14:51:30 |
| pass3 | DOAJ | W-type strontium hexaferrite dielectric | 0 | 14:51:30 |
| pass3 | SemanticScholar | W-type strontium hexaferrite dielectric | 0 | 14:51:30 |
| pass3 | EuropePMC | Mossbauer Debye temperature recoil-free fraction ferrite sites | 4 | 14:51:30 |
| pass3 | arXiv | Mossbauer Debye temperature recoil-free fraction ferrite sites | 0 | 14:51:30 |
| pass3 | Crossref-query | Mossbauer Debye temperature recoil-free fraction ferrite sites | 1747830 | 14:51:30 |
| pass3 | SemanticScholar | Mossbauer Debye temperature recoil-free fraction ferrite sites | 0 | 14:51:30 |
| pass3 | EuropePMC | recoil-free fraction Fe sites hexaferrite Mossbauer | 1 | 14:51:30 |
| pass3 | arXiv | recoil-free fraction Fe sites hexaferrite Mossbauer | 0 | 14:51:30 |
| pass3 | Crossref-query | recoil-free fraction Fe sites hexaferrite Mossbauer | 1124052 | 14:51:30 |
| pass3 | SemanticScholar | recoil-free fraction Fe sites hexaferrite Mossbauer | 0 | 14:51:30 |
| pass3 | EuropePMC | Schlomann porous ferrite linewidth | 0 | 14:51:30 |
| pass3 | arXiv | Schlomann porous ferrite linewidth | 0 | 14:51:30 |
| pass3 | Crossref-query | Schlomann porous ferrite linewidth | 240650 | 14:51:30 |
| pass3 | SemanticScholar | Schlomann porous ferrite linewidth | 0 | 14:51:30 |
| pass3 | EuropePMC | ferromagnetic resonance hexaferrite Kittel uniaxial anisotropy | 5 | 14:51:30 |
| pass3 | arXiv | ferromagnetic resonance hexaferrite Kittel uniaxial anisotropy | 0 | 14:51:30 |
| pass3 | Crossref-query | ferromagnetic resonance hexaferrite Kittel uniaxial anisotropy | 484226 | 14:51:30 |
| pass3 | SemanticScholar | ferromagnetic resonance hexaferrite Kittel uniaxial anisotropy | 0 | 14:51:30 |
| pass3 | EuropePMC | Smit Beljers resonance condition | 8 | 14:51:30 |
| pass3 | arXiv | Smit Beljers resonance condition | 0 | 14:51:30 |
| pass3 | Crossref-query | Smit Beljers resonance condition | 548780 | 14:51:30 |
| pass3 | SemanticScholar | Smit Beljers resonance condition | 0 | 14:51:30 |
| pass3 | EuropePMC | ferromagnetic resonance polycrystalline hexagonal ferrite anisotropy field distribution | 46 | 14:51:30 |
| pass3 | arXiv | ferromagnetic resonance polycrystalline hexagonal ferrite anisotropy field distribution | 0 | 14:51:30 |
| pass3 | Crossref-query | ferromagnetic resonance polycrystalline hexagonal ferrite anisotropy field distribution | 2059423 | 14:51:30 |
| pass3 | SemanticScholar | ferromagnetic resonance polycrystalline hexagonal ferrite anisotropy field distribution | 0 | 14:51:30 |
| pass3 | EuropePMC | Mossbauer spectra ferrimagnetic hexagonal ferrite sublattice assignment | 2 | 14:51:30 |
| pass3 | arXiv | Mossbauer spectra ferrimagnetic hexagonal ferrite sublattice assignment | 0 | 14:51:30 |
| pass3 | Crossref-query | Mossbauer spectra ferrimagnetic hexagonal ferrite sublattice assignment | 276984 | 14:51:30 |
| pass3 | SemanticScholar | Mossbauer spectra ferrimagnetic hexagonal ferrite sublattice assignment | 0 | 14:51:30 |
| pass4 | arXiv | all:hexaferrite | 112 | 15:11:31 |
| pass4 | arXiv | all:"W-type" AND all:ferrite | 2 | 15:11:31 |
| pass4 | arXiv | all:"W-type" AND all:hexaferrite | 4 | 15:11:31 |
| pass4 | arXiv | all:hexagonal AND all:ferrite AND all:anisotropy | 6 | 15:11:31 |
| pass4 | arXiv | all:SrFe12O19 | 12 | 15:11:31 |
| pass4 | arXiv | all:BaFe12O19 | 24 | 15:11:31 |
| pass4 | arXiv | all:"M-type" AND all:hexaferrite AND all:DFT | 0 | 15:11:31 |
| pass4 | arXiv | all:Mossbauer AND all:hexaferrite | 1 | 15:11:31 |
| pass4 | arXiv | all:"ferromagnetic resonance" AND all:hexaferrite | 5 | 15:11:31 |
| pass4 | arXiv | all:"strontium hexaferrite" | 15 | 15:11:31 |
| pass4 | arXiv | all:"Fe18O27" | 0 | 15:11:31 |
| pass4 | arXiv | all:"first-principles" AND all:hexaferrite | 7 | 15:11:31 |
| pass4 | arXiv | all:"recoil-free" AND all:Mossbauer AND all:ferrite | 0 | 15:11:31 |

Пояснение: pass=проход 1 (OpenAlex+Crossref bibliographic), pass2 (Europe PMC, arXiv, Crossref query, DOAJ), pass3 (то же + Semantic Scholar, DFT/методология), pass4 (arXiv, простые запросы, до 120 записей на запрос). Нулевые значения у OpenAlex в pass1 после 21-го запроса и у Semantic Scholar -- отказ по лимиту.

## Статус L04 / L06 / L07 / L08 (OA-версия)
| ID | DOI | OA найдена? | чем проверено |
|---|---|---|---|
| L04 Islam 2024 PRB | 10.1103/PhysRevB.109.024414 | нет | Semantic Scholar openAccessPdf=CLOSED; Crossref лицензия APS default; arXiv (all:SrFe18O27, au:Christensen, W-type) -- 0 записей; OpenAlex is_oa=False (проход 1) |
| L06 Zn2-xCoxW FMR 1995 | 10.1016/0304-8853(95)00406-8 | нет | S2 CLOSED; Elsevier TDM-лицензия (не OA); OpenAlex is_oa=False; arXiv нет |
| L07 Rancourt 1989 | 10.1016/0168-583X(89)90428-X | нет | S2 CLOSED; Elsevier TDM; Europe PMC: нет |
| L08 Newbury & Ritchie 2013 | 10.1002/sca.21041 | нет | S2 CLOSED; Wiley VOR (не OA); Europe PMC: не в PMC |

## Связи с корпусом
В candidates.csv колонка corpus_link: найдены L01 (C017), L02 (C009), L03 (C010), L04 (C006), L05 (C003), L06 (C068), L11 (C090), плюс L07 (C134) и L08 (C135), добавленные по DOI. Не найдены поисковыми запросами (не гексаферритная тема или нет DOI): L09, L10, L12, L13, L14.

## Воспроизведение
Скрипты в scripts\: search.py, search2.py, search3.py, arxiv_sweep.py (поиск), merge.py, build.py (отбор и relevance), oa.py, download.py, download2.py, epmc_xml.py (OA и загрузка), final.py (CSV). Правила relevance: 1 = Sr W-тип (по заголовку/реферату); 2 = прочие W-тип с свойствами/Co,Zn,Mn/Ba, плюс исторические работы по анизотропии/Мёссбауэру/нейтронам; 3 = W-тип прочее (DFT W, диэлектрика, прикладные); 4 = методологические аналоги (M-тип DFT, FMR, Мёссбауэр). Оценка эвристическая (регулярные выражения по заголовку и реферату) с ручными правками; не заменяет чтения. 11 записей без реферата и с общим названием отсеяно в results\screened_out.csv.
