import React, { useState } from 'react';
import { styled } from '@mui/material/styles';
import { Link as NavLink, useRouteMatch } from 'react-router-dom';
import Box from '@mui/material/Box';
import MuiDrawer from '@mui/material/Drawer';
import Stack from '@mui/material/Stack';
import List from '@mui/material/List';
import CssBaseline from '@mui/material/CssBaseline';
import AccountBalanceIcon from '@mui/icons-material/AccountBalance';
import ListItemButton from '@mui/material/ListItemButton';
import ListItemIcon from '@mui/material/ListItemIcon';
import ListItemText from '@mui/material/ListItemText';
import LiveHelpIcon from '@mui/icons-material/LiveHelp';
import Avatar from '@mui/material/Avatar';
import MenuBookIcon from '@mui/icons-material/MenuBook';
import MapIcon from '@mui/icons-material/Map';
import TaskIcon from '@mui/icons-material/Task';
import { Typography } from '@mui/material';
import profiledefault from '../../../assets/user.png';
import logo from '../../../assets/logo.svg';
import { SettingsDropdown } from '../NavigationBar/SettingsDropdown/SettingsDropdown';
import styles from './drawer.module.css';

const drawerWidth = 240;

const openedMixin = (theme) => ({
  width: drawerWidth,
  transition: theme
    .transitions
    .create('width', {
      easing: theme.transitions.easing.smooth,
      duration: theme.transitions.duration.enteringScreen,
    }),
  overflowX: 'hidden',
});

const closedMixin = (theme) => ({
  transition: theme
    .transitions
    .create('width', {
      easing: theme.transitions.easing.smooth,
      duration: theme.transitions.duration.leavingScreen,
    }),
  overflowX: 'hidden',
  [
  theme
    .breakpoints
    .up('sm')
  ]: {
    width: `calc(${theme.spacing(10)} + 1px)`,
  },
});

const DrawerHeader = styled('div')(({ theme }) => ({
  display: 'flex',
  alignItems: 'center',
  justifyContent: 'flex-end',
  padding: theme.spacing(0, 1),
}));

const Drawer = styled(MuiDrawer, {
  shouldForwardProp: (prop) => prop !== 'open',
})(({ theme, open }) => ({
  width: drawerWidth,
  flexShrink: 0,
  whiteSpace: 'nowrap',
  boxSizing: 'border-box',
  ...(open && {
    ...openedMixin(theme),
    '& .MuiDrawer-paper': openedMixin(theme),
  }),
  ...(!open && {
    ...closedMixin(theme),
    '& .MuiDrawer-paper': closedMixin(theme),
  }),
}));

function indexIcon(index) {
  switch (index) {
    case 0:
      return (
        <TaskIcon
          sx={{
            color: '#00579b',
          }}
        />
      );
    case 1:
      return (
        <AccountBalanceIcon
          sx={{
            color: '#00579b',
          }}
        />
      );
    case 2:
      return (
        <MapIcon
          sx={{
            color: '#00579b',
          }}
        />
      );
    case 3:
      return (
        <MenuBookIcon
          sx={{
            color: '#00579b',
          }}
        />
      );
    default:
      return (
        <LiveHelpIcon
          sx={{
            color: '#00579b',
          }}
        />
      );
  }
}

export default function MiniDrawer({ user, setUser }) {
  const [open,
    setOpen] = useState(false);
  const [dropDown, setDropDown] = useState(false);

  const toggleDrawer = (isOpen) => {
    setOpen(isOpen);
    if (!isOpen) setDropDown(false);
  };

  const onMouseEnter = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(true);
  };

  const { url } = useRouteMatch();
  // TODO: change routes after implementation
  const routes = ['dashboard', 'dashboard', 'dashboard', 'documentation', 'help'];
  const titles = ['Projects', 'Institutions', 'Gene Mapper', 'Documentation', 'Help'];

  return (
    <Box sx={{
      display: 'flex',
    }}
    >
      <CssBaseline />
      <Drawer
        PaperProps={{
          className: styles.drawer,
        }}
        variant="permanent"
        open={open}
        onMouseOver={() => toggleDrawer(true)}
        onMouseLeave={() => toggleDrawer(false)}
      >
        <DrawerHeader>
          <div
            className={styles.drawerHeader}
          >
            <Avatar
              alt="logo"
              src={logo}
              className={styles.avatar}
            />
          </div>
        </DrawerHeader>
        <List className={styles.drawerList}>
          {titles.map((text, index) => (
            <NavLink
              to={`${url}/${routes[index]}`}
              className={styles.drawerNavLink}
            >
              <ListItemButton
                key={text}
                onMouseOver={() => setDropDown(false)}
                sx={{
                  minHeight: 48,
                  padding: '10px',
                  '&:hover': {
                    backgroundColor: '#01579B',
                  },
                }}
              >
                <div className={styles.navbarItemContainer}>
                  <ListItemIcon
                    className={styles.drawerListItemIcon}
                  >
                    {indexIcon(index)}
                  </ListItemIcon>
                </div>
                <ListItemText
                  primary={text}
                  sx={{
                    opacity: open
                      ? 1
                      : 0,
                    px: 1,
                  }}
                />
              </ListItemButton>
            </NavLink>
          ))}
        </List>
        <Stack
          onMouseEnter={onMouseEnter}
          className={styles.profileAvatarContainer}
        >
          {dropDown
          && <SettingsDropdown setUser={setUser} handleClose={() => setDropDown(false)} />}
          <Avatar
            alt="profile"
            src={profiledefault}
            className={styles.profileAvatar}
          />
          <Typography
            sx={{
              opacity: open
                ? 1
                : 0,
              color: 'white',
            }}
          >
            {`Hi, ${user.firstName}!`}
          </Typography>
        </Stack>
      </Drawer>
    </Box>
  );
}
