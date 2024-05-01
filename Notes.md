#### Open ports
```
rustscan -a $ip --ulimit 5000
```

#### IIS Webdav
```
cadaver $ip
upload conptyshell
```

#### LDAP
```
nxc ldap $ip -u '' -p '' -M get-desc-user
```

#### enum4linux-ng
```
enum4linux-ng -A $ip
```

#### nxc
```
nxc smb $ip -u 'guest' -p '' --spider C\$ --pattern txt
nxc smb $ip -u 'guest' -p '' -M spider_plus -o DOWNLOAD_FLAG=True
nxc smb $ip -u 'guest' -p '' --users
```