import click
from classes import *
from fgraph import *

@click.group()
def cli():
    """
    Основная команда программы. Используйте команды ниже для управления программой.
    """
    pass


@cli.command()
@click.argument('model_b', type=click.Choice(['bgi', 'mhd']))
@click.argument('model_psi', type=click.Choice(['old', 'new']))
@click.option('--sargs', '-sa', nargs=2, default=(1, 1), type=click.Tuple([int, int]), help='Дополнительные аргументы системы')
@click.option('--system', '-sys', is_flag=True, help='Просмотр параметров системы')
@click.option('--allpulsar', '-all', is_flag=True, help='Просмотр параметров всех пульсаров')
@click.option('--npulsar', '-puls', default=1, type=int, help='Просмотр параметров пульсара')
@click.option('--keysys', '-ks', multiple=True, default=('model_b','model_psi', 'psi_args'), type=click.Choice(['model_b', 'model_psi', 'num_pulsars', 
                                                                                                                'num_points', 'psi_args', 'b_args', 
                                                                                                                'other_args']), help='Параметры для пульсара')
@click.option('--keypuls', '-kp', multiple=True, default=('P', 'P_dot', 'chi_d'), type=click.Choice(['name', 'chi','chi_d', 'beta', 
                                                                                                     'beta_d','P','P_dot', 'Omega', 'model_args', 
                                                                                                     'W50','sigma', 'R_0', 'B', 'phi_m']), help='Параметры для пульсара')
def check(model_b, model_psi, sargs, system, allpulsar, npulsar, keypuls, keysys):
    """
    Проверка параметров системы и пульсара.

    :param model_b: модель магнитного поля (bgi или mhd)

    :param model_psi: модель поверхностной магнитной полярности (old или new)
    """
    system_p = System(model_b, model_psi, 1, sargs)
    table = Table(header_style="cyan")
    console = Console()
    
    if system:
        table.add_column("parameters", justify = "centr", style="cyan")
        table.add_column("value", justify = "center", style="cyan")
        for key in keysys:
            table.add_row(key, str(system_p.const[key]))
        console.print('')
        console.print(table)
    elif allpulsar:
        table.add_column("pulsar", justify = "centr", style="cyan")
        for key in keypuls:
            table.add_column(key, justify = "centr", style="cyan")
        for k in range(0,8):
            values = tuple([str(system_p.pulsars[k].const[key]) for key in keypuls])
            table.add_row(system_p.pulsars[k].name, *values)
        console.print('')
        console.print(table)
    else:
        table.add_column("parameter", justify = "centr", style="cyan")
        for key in keypuls:
            table.add_column(key, justify = "centr", style="cyan")
        values = tuple([str(system_p.pulsars[npulsar].const[key]) for key in keypuls])
        table.add_row(system_p.pulsars[npulsar].name, *values)
        console.print('')
        console.print(table)


@cli.command()
@click.argument('parameter', type=click.Choice(['gamma', 'hgap', 'psi', 'field', 'lambda', 'prof']))
@click.argument('form', type=click.Choice(['save', 'print', 'table']))
@click.option('--bmodel', '-b', default='mhd', type=click.Choice(['bgi', 'mhd']), help='Модель магнитного поля')
@click.option('--psimodel', '-psi', default='old', type=click.Choice(['old', 'new']), help='Модель электрического поля')
@click.option('--npoints', '-num', default=100, type=int, help="Количество точек на графике")
@click.option('--npulsar', '-puls', default=1, type=int, help='Номер пульсара')
@click.option('--allpulsar', '-all', is_flag=True, help='Все пульсары')
@click.option('--sargs', '-sa', nargs=2, default=(5, 5), type=click.Tuple([int, int]), help='Дополнительные аргументы системы')
@click.option('--args', '-a', nargs=2, default=(0.5, 45), type=click.Tuple([float, int]), help='Дополнительные аргументы параметра')
def test(parameter, form, bmodel, psimodel, npoints, npulsar, allpulsar, sargs, args):
    """
    Изучение физических зависимостей.

    :param parameter: параметр для тестирования
    :param form: формат вывода результатов
    """
    system = System(bmodel, psimodel, npoints, sargs)
    
    if allpulsar:
        num_pulsars = system.const['num_pulsars']
        for k_pulsar in range(num_pulsars):
            pulsar_num = k_pulsar + 1
            test_parameter(system, pulsar_num, parameter, form, args)
    else:
        pulsar_num = npulsar
        test_parameter(system, pulsar_num, parameter, form, args)


def test_parameter(system, pulsar_num, parameter, form, args):
    """
    Выполняет тестирование параметра для указанного пульсара.

    :param system: система пульсара
    :param pulsar_num: номер пульсара
    :param parameter: параметр для тестирования (gamma, hgap, psi, field, lambda, prof)
    :param form: формат вывода результатов (table или plot)
    :param args: дополнительные аргументы
    """
    if parameter == "gamma":
        r, phi = args
        H, Gamma, Delta = system.get_gamma(pulsar_num - 1, float(r), int(phi))
        graph_gamma(H, Gamma, Delta, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b, 'r': str(r), 'phi': str(phi)})
    elif parameter == "field":
        r, phi = args
        H, Field = system.get_field(pulsar_num - 1, float(r), float(phi))
        graph_field(H, Field, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b, 'r': str(r), 'phi': str(phi)})
    elif parameter == "hgap":
        H, Hgap = system.get_hgap(pulsar_num - 1)
        graph_hgap(H, Hgap, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b})
    elif parameter == "psi":
        R, Psi = system.get_psi(pulsar_num - 1)
        graph_psi(R, Psi, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b})
    elif parameter == "lambda":
        X, Y, Lambda = system.get_lambda(pulsar_num - 1)
        graph_lambda(X, Y, Lambda, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b})
    elif parameter == "prof":
        R, LambdaMP, LambdaIP = system.get_profile(pulsar_num - 1)
        graph_profile(R, LambdaMP, LambdaIP, form, {'name': system.pulsars[pulsar_num - 1].name, 'model_b': system.model_b})
    else:
        click.echo("CommandError")


if __name__ == '__main__':
    cli()