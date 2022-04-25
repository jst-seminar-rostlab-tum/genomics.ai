import React from 'react';
import LiveHelpIcon from '@mui/icons-material/LiveHelp';
import MenuBookIcon from '@mui/icons-material/MenuBook';
import MapIcon from '@mui/icons-material/Map';
import TaskIcon from '@mui/icons-material/Task';
import AccountBalanceIcon from '@mui/icons-material/AccountBalance';
import List from '@mui/material/List';
import Tooltip from '@mui/material/Tooltip';
import Box from '@mui/material/Box';
import { Link as NavLink, useRouteMatch, useLocation } from 'react-router-dom';
import ListItemIcon from '@mui/material/ListItemIcon';
import LogoutIcon from '@mui/icons-material/Logout';
import SettingsIcon from '@mui/icons-material/Settings';
import geneIcon from 'assets/gene.png';
import styles from './sidebar.module.css';

function indexIcon(index) {
  switch (index) {
    case 0:
      return (
        <img
          alt="gene-icon"
          src={geneIcon}
          className={styles.geneIcon}
        />
      );
    case 1:
      return (<TaskIcon className={styles.coloredIcon} />);
    case 2:
      return (<AccountBalanceIcon className={styles.coloredIcon} />);
    case 3:
      return (<MapIcon className={styles.coloredIcon} />);
    case 4:
      return (<MenuBookIcon className={styles.coloredIcon} />);
    default:
      return (<LiveHelpIcon className={styles.coloredIcon} />);
  }
}

export default function Sidebar(props) {
  const { setUser } = props;
  const routes = ['dashboard', 'teams', 'institutions', 'genemapper', 'documentation', 'help'];
  const titles = ['Dashboard', 'Teams', 'Institutions', 'Gene Mapper', 'Documentation', 'Help'];
  const { url } = useRouteMatch();
  const location = useLocation();
  const path = location.pathname;
  const settingsPath = '/sequencer/settings';

  return (
    <Box>
      <Box className={styles.sidebarNav}>
        <Box className={styles.sidebarWrap}>
          <List className={styles.iconList}>
            {routes.map((route, index) => (
              <NavLink
                className={styles.navlink}
                to={`${url}/${route}`}
                key={route.toString()}
              >
                <Tooltip
                  title={titles[index]}
                  placement="right"
                  componentsProps={{
                    tooltip: {
                      sx: {
                        bgcolor: '#5676E4',
                      },
                    },
                  }}
                >
                  <Box
                    className={styles.navbarItemContainer}
                    sx={{ background: path.includes(route) ? '#5676E5' : '#184060' }}
                  >
                    <ListItemIcon className={styles.listItemIcon}>
                      {indexIcon(index)}
                    </ListItemIcon>
                  </Box>
                </Tooltip>
              </NavLink>
            ))}
            <NavLink
              to={settingsPath}
              className={`${styles.navlinkIcon} ${styles.bottomIcons} ${styles.settingsIcon}`}
            >
              <SettingsIcon />
            </NavLink>
            <NavLink
              to="/"
              onClick={() => {
                setUser(null);
                localStorage.removeItem('user');
                localStorage.removeItem('jwt');
              }}
              className={`${styles.navlinkIcon} ${styles.bottomIcons}`}
            >
              <LogoutIcon />
            </NavLink>
          </List>
        </Box>
      </Box>
    </Box>
  );
}
