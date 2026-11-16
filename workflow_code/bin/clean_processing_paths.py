#!/usr/bin/env python3
"""Remove local directory prefixes from a copy of the processing records."""

import argparse
import json
import re
from pathlib import Path


def path_pattern(prefixes):
    prefixes = {str(prefix).rstrip('/') for prefix in prefixes}
    if not prefixes:
        return None
    if any(not prefix.startswith('/') for prefix in prefixes):
        raise ValueError('Path prefixes must be absolute directories other than /')
    original_prefixes = list(prefixes)
    prefixes.update(prefix.replace("'", "'\"'\"'") for prefix in original_prefixes)
    prefixes.update(re.sub(r'''([*?{}\[\]'" ()\\!&|<>`:])''', r'\\\1', prefix)
                    for prefix in original_prefixes)
    prefixes = sorted(prefixes, key=len, reverse=True)
    return re.compile(r'''(?<![^\s="'(])(?:''' + '|'.join(re.escape(prefix) for prefix in prefixes)
                      + r''')(?P<separator>/|(?=$|[\s"'<>|;(),\\]))''')


def clean(text, pattern):
    if pattern is None:
        return text
    return pattern.sub(lambda match: '.' + match['separator'], text)


def clean_json(value, pattern):
    if isinstance(value, str):
        return clean(value, pattern)
    if isinstance(value, list):
        return [clean_json(item, pattern) for item in value]
    if isinstance(value, dict):
        return {key: clean_json(item, pattern) for key, item in value.items()}
    return value


def recorded_prefixes(directory):
    prefixes = []
    for path in directory.rglob('workflow_run.json'):
        record = json.loads(path.read_text())
        for key in ('launch_dir', 'project_dir', 'work_dir', 'output_dir'):
            value = record.get(key)
            if isinstance(value, str) and value.startswith('/') and value != '/':
                prefixes.append(value)
    return prefixes


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('directory', type=Path)
    parser.add_argument('--prefix', action='append', default=[], help='Local directory prefix to remove (repeatable)')
    args = parser.parse_args()
    try:
        pattern = path_pattern(args.prefix + recorded_prefixes(args.directory))
    except ValueError as error:
        parser.error(str(error))
    for path in args.directory.rglob('*'):
        if path.is_file():
            text = path.read_text()
            if path.suffix == '.json':
                text = json.dumps(clean_json(json.loads(text), pattern), indent=4) + '\n'
            else:
                text = clean(text, pattern)
            path.write_text(text)


if __name__ == '__main__':
    main()
