package Cache;
use strict; use warnings; use 5.008;
use Digest::SHA qw(sha1_hex);
use File::Path qw(mkpath);
use File::Spec;

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

sub try_get {
    my $self = shift;
    my $hash = sha1_hex(@_);
    $self->{cache_file} =
      File::Spec->catfile($temp_dir, $self->{cache_dir}, $hash);

    my $file = $self->{cache_file};
    return '' unless -e $file;

    # Oh look, it exists.
    {
        local $/;
        open(my $h, $file)
          or die "Cache: unable to open cache \"$file\"";
        return <$h>;
    }
}

sub new_store {
    my $self = shift;
    my $str = shift;
    my $file = $self->{cache_file};

    # We'll just silently fail for now.
    open(my $h, ">$file") or return;
    print $h $str;
    close($h);
}

1;
