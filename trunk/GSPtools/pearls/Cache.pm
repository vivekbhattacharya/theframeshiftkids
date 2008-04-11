package Cache;
use strict; use warnings;
use Digest::SHA qw(sha1_hex);
use File::Spec;
use File::Path qw(mkpath);

my $temp_dir = File::Spec->tmpdir();

sub init {
    my $self = {};
    my ($class, $dir) = @_;
    $self->{cache_dir} = $dir;
    bless($self, $class);

    $self->make_cache_dir;
    return $self;
}

sub make_cache_dir {
    my $self = shift;
    my $dir = File::Spec->catdir($temp_dir, $self->{cache_dir});
    return if -d $dir;
    mkpath $dir;
}

sub maybe_print_and_exit {
    my $self = shift;
    my $hash = sha1_hex(@_);
    $self->{cache_file} =
      File::Spec->catfile($temp_dir, $self->{cache_dir}, $hash);

    my $file = $self->{cache_file};
    return if !-e $file;

    # Oh look, it exists, we can quit.
    open(my $h, $file);
    print while <$h>;
    exit;
}

sub new_store {
    my $self = shift;
    my $str = shift;
    my $file = $self->{cache_file};
    open(my $h, ">$file")
      or die 'Cache.pm: Unable to write to store.';
    print $h $str;
    close($h);
}

1;
