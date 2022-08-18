import configargparse
import argparse


class write_config(object):
    def __init__(self, namespace):
        self.namespace = namespace

    def write_action(self, a):

        is_ = False
        for item in a.option_strings:
            is_ = is_ or item.startswith('--')
            is_ = is_ and item != '--help'
        is_ = is_ and not getattr(a, "not_in_write_out_config_file", False)
        is_ = is_ and not getattr(a, "is_config_file_arg", False)
        #is_ = is_ and not getattr(a, "is_write_out_config_file_arg", False)

        if is_:
            msg = '#help: '
            msg += f"{', '.join(a.option_strings)}; "
            msg += f"default: {a.default}; " if a.default is not None else 'default: ;'
            if a.choices is not None:
                msg += f"choices: {', '.join(map(str,a.choices))}; "
            msg += '\n# '.join(a.help.split('\n'))
            value = getattr(self.namespace, a.dest, '')
            if isinstance(value, list):
                value = "["+", ".join(map(str, value))+"]"
            if value is None:
                msg += f"\n#{a.dest} = "
            elif isinstance(a.const, bool):
                msg += f"\n{a.dest}" if a.const == value else f"\n#{a.dest}"
            elif isinstance(value, bool) and a.default.lower() in ['yes', 'no']:
                msg += f"\n{a.dest} = yes" if value else f"\n{a.dest} = no"
            else:
                msg += f"\n{a.dest} = {value}"
            return msg+'\n'
        else:
            return ''

    def write_group(self, g):
        msg = ''
        for a in g._group_actions:
            msg += self.write_action(a)
        if msg != '':
            title = '\n## '.join(g.title.split('\n'))
            msg = f"## {title}\n" + msg + '\n'
        return msg

    def write_groups(self, gs):
        msg = ''
        for g in gs:
            msg += self.write_group(g)
        return msg

def add_argument(self, *args, **kwargs):
    """
    Arguments added:
        not_in_write_out_config_file: excluded in 'is_write_out_config_file_arg' file
        not_add_hyphen_option: bool
    """
    not_in_write_out_config_file = kwargs.pop("not_in_write_out_config_file", None)
    not_add_hyphen_option = kwargs.pop("not_add_hyphen_option", False)
    
    if not not_add_hyphen_option:
        args = list(args)
        for op in args:
            if '_' in op and op.startswith('-'):
                op_hy = op.replace('_', '-')
                if op_hy not in args:
                    args += [op_hy, ]
    
    action = self.configargparse_original_add_argument(*args, **kwargs)
    action.not_in_write_out_config_file = not_in_write_out_config_file
    
    
    
    return action


add_argument.__doc__ += configargparse.ArgParser.add_argument.__doc__


class ArgumentParser(configargparse.ArgParser):
    def write_config_file(self, parsed_namespace, output_file_paths, exit_after=False):
        if output_file_paths:
            ww = write_config(parsed_namespace)
            for output_file_path in output_file_paths:
                with open(output_file_path,'w') as f:
                    f.write(ww.write_groups(self._action_groups))
            print("Wrote config file to " + ", ".join(output_file_paths))
            if exit_after:
                self.exit(0)
                
        # del argument with not_in_write_out_config_file in parsed_namespace
        #[delattr(parsed_namespace, a.dest) for a in self._actions if getattr(a, "not_in_write_out_config_file", False)]
        #super().write_config_file(parsed_namespace, output_file_paths, exit_after)
    #write_config_file.__doc__ = configargparse.ArgParser.write_config_file.__doc__

argparse._ActionsContainer.configargparse_original_add_argument = configargparse.ArgParser.add_argument
argparse._ActionsContainer.add_argument = add_argument
argparse._ActionsContainer.add_arg = add_argument
argparse._ActionsContainer.add = add_argument

ArgParser = ArgumentParser
Parser = ArgumentParser
