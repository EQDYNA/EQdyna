#! /usr/bin/env python3
"""
Fetch EVERY public SCEC/USGS cvws submission belonging to the EQdyna owner
(`dliu`, `dliu.2`, ...) across ALL benchmarks, and lay them out in
scec_archive/<benchmark>/<code-version>-<resolution>-<year>/ per
scec_archive/README.md.

Request shape is the one proven by scratch/tpv29/scoring/fetch_cvws_tpv29.py,
learned by scraping the CGI's own forms:

  step 1  benchmark list  -F o=1005 -F 'G0012=Go -->'
  step 2  user list       -F o=1005 -F 'G1045<bm>= Select '   -> hidden Q0001
  step 3  file list       -F o=1005 -F m=<bm> -F Q0001=<users>
                          -F 'G1047<user>= Select '
                          -> hidden Q0001 = on-fault stations (*-joined)
                             hidden Q0003 = off-fault stations (*-joined)
                             one 'Raw Data' submit button per served file
  step 4  raw file        -F o=1005 -F m=<bm> -F mus=<user>
                          -F Q0001=<on-fault> -F Q0003=<off-fault>
                          -F '<button name>=Raw Data'

The submit-button NAME is per-file and is read off the step-3 page, never
guessed:  G1063cplot for cplot, G1059<station> for faultst*/body*.

The <pre> rendering prepends ONE space to lines starting with '#'; unrender()
reverses it, giving the file exactly as uploaded.
"""
import argparse, datetime, html, json, os, re, subprocess, sys, time

CGI = 'https://strike.scec.org/cvws/cgi-bin/cvws.cgi'
HERE = os.path.dirname(os.path.abspath(__file__))
ARCHIVE = os.path.abspath(os.path.join(HERE, '..', '..', 'scec_archive'))
TODAY = datetime.date.today().isoformat()


def curl(args, timeout=300, tries=3):
    last = None
    for i in range(tries):
        r = subprocess.run(['curl', '-sS', '-m', str(timeout)] + args,
                           capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout
        last = r.stderr.strip()
        time.sleep(2 * (i + 1))
    raise RuntimeError(f'curl failed after {tries}: {last}')


def scrape_forms(doc):
    forms = []
    for f in re.findall(r'(?is)<form[^>]*>(.*?)</form>', doc):
        fields = []
        for tag, attrs in re.findall(
                r'(?is)<(input|select|textarea|button)\b([^>]*)>', f):
            a = {k.lower(): html.unescape(v)
                 for k, v in re.findall(r'(?is)(\w+)\s*=\s*"([^"]*)"', attrs)}
            a['tag'] = tag.lower()
            fields.append(a)
        forms.append(fields)
    return forms


def hidden(forms, name):
    for ff in forms:
        for a in ff:
            if a.get('type') == 'hidden' and a.get('name') == name:
                return a.get('value', '')
    return None


def benchmarks():
    doc = curl(['-F', 'o=1005', '-F', 'G0012=Go -->', CGI])
    return sorted(set(re.findall(r'name="G1045([^"]+)"', doc)))


def user_page(bm):
    return curl(['-F', 'o=1005', '-F', f'G1045{bm}= Select ', CGI])


def users_and_labels(doc):
    """(-> list of users, {user: descriptive label}) from a step-2 page."""
    ul = hidden(scrape_forms(doc), 'Q0001')
    users = [u for u in (ul or '').split('*') if u]
    labels = {}
    txt = html.unescape(re.sub(r'(?is)<t[dh][^>]*>', '\x01', doc))
    txt = re.sub(r'<[^>]+>', ' ', txt)
    cells = [re.sub(r'\s+', ' ', c).strip() for c in txt.split('\x01')]
    for i, c in enumerate(cells):
        if c in users and i + 1 < len(cells):
            labels.setdefault(c, cells[i + 1])
    return users, labels, (ul or '')


def discover(bm, ulist, user, dump=None):
    doc = curl(['-F', 'o=1005', '-F', f'm={bm}', '-F', f'Q0001={ulist}',
                '-F', f'G1047{user}= Select ', CGI])
    forms = scrape_forms(doc)
    if dump:
        json.dump({'benchmark': bm, 'user': user, 'url': CGI, 'forms': forms},
                  open(dump, 'w'), indent=1)
    btn = {}
    for ff in forms:
        for a in ff:
            n, v = a.get('name', ''), a.get('value', '').strip()
            if v == 'Raw Data' and n:
                # strip the known G#### prefix -> the file name
                m = re.match(r'^(G\d{4})(.+)$', n)
                if m:
                    btn[m.group(2)] = n
    return dict(buttons=btn,
                q0001=hidden(forms, 'Q0001') or '',
                q0003=hidden(forms, 'Q0003') or '')


def unrender(text):
    return '\n'.join(l[1:] if l.startswith(' #') else l
                     for l in text.split('\n'))


def fetch_file(bm, user, field, q1, q3):
    doc = curl(['-F', 'o=1005', '-F', f'm={bm}', '-F', f'mus={user}',
                '-F', f'Q0001={q1}', '-F', f'Q0003={q3}',
                '-F', f'{field}=Raw Data', CGI])
    pre = re.findall(r'(?is)<pre[^>]*>(.*?)</pre>', doc)
    if not pre:
        why = ('Data File Not Found' if 'Data File Not Found' in doc
               else f'no <pre>, {len(doc)} B')
        raise RuntimeError(why)
    return unrender(html.unescape(re.sub(r'<[^>]+>', '', pre[0])))


HDR = re.compile(r'^#\s*(\w+)\s*=\s*(.*?)\s*$')


def headers(text):
    out = {}
    for line in text.split('\n')[:40]:
        if not line.startswith('#'):
            if out:
                break
            continue
        m = HDR.match(line)
        if m:
            out.setdefault(m.group(1), m.group(2))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--who', default='dliu')
    ap.add_argument('--only', nargs='*', help='limit to these benchmarks')
    ap.add_argument('--enumerate-only', action='store_true')
    ap.add_argument('--out', default=os.path.join(HERE, 'manifest.json'))
    a = ap.parse_args()

    bms = a.only or benchmarks()
    print(f'{len(bms)} benchmarks offered by the public area', file=sys.stderr)

    found = {}
    for bm in bms:
        try:
            doc = user_page(bm)
            users, labels, ulist = users_and_labels(doc)
        except Exception as e:
            print(f'  {bm}: {e}', file=sys.stderr)
            continue
        mine = [u for u in users
                if u == a.who or u.startswith(a.who + '.')]
        if mine:
            found[bm] = dict(users=mine,
                             labels={u: labels.get(u, '') for u in mine},
                             ulist=ulist, all_users=users)
            print(f'  {bm}: {mine}', file=sys.stderr)
    json.dump(found, open(os.path.join(HERE, 'enumeration.json'), 'w'), indent=1)
    if a.enumerate_only:
        print(json.dumps({k: v['labels'] for k, v in found.items()}, indent=1))
        return

    manifest = []
    for bm, info in sorted(found.items()):
        for user in info['users']:
            d = discover(bm, info['ulist'], user,
                         dump=os.path.join(HERE, 'forms',
                                           f'{bm}_{user}.json'))
            files = sorted(d['buttons'])
            print(f'== {bm}/{user}: {len(files)} files offered',
                  file=sys.stderr)
            got, failed = {}, {}
            texts = {}
            for f in files:
                try:
                    t = fetch_file(bm, user, d['buttons'][f],
                                   d['q0001'], d['q0003'])
                except Exception as e:
                    failed[f] = str(e)
                    print(f'   !! {f}: {e}', file=sys.stderr)
                    continue
                texts[f] = t
                got[f] = len(t)
            manifest.append(dict(benchmark=bm, user=user,
                                 label=info['labels'].get(user, ''),
                                 buttons=d['buttons'], got=got, failed=failed,
                                 headers=headers(texts.get('cplot',
                                                 next(iter(texts.values()), ''))),
                                 fetched=TODAY))
            # stash raw texts for the layout pass
            sd = os.path.join(HERE, 'raw', bm, user)
            os.makedirs(sd, exist_ok=True)
            for f, t in texts.items():
                open(os.path.join(sd, f), 'w').write(t)
    json.dump(manifest, open(a.out, 'w'), indent=1)
    print(f'wrote {a.out}', file=sys.stderr)


if __name__ == '__main__':
    os.makedirs(os.path.join(HERE, 'forms'), exist_ok=True)
    main()
