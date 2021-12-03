import React, { useState } from 'react';
import { Link as NavLink, useRouteMatch } from 'react-router-dom';
import { SettingsDropdown } from './SettingsDropdown/SettingsDropdown';
import styles from './navigationBar.module.css';
import geneIcon from '../../../assets/logo-blue.png';
import profiledefault from '../../../assets/user.png';

function NavigationBar({ sidebarShown, user, setUser }) {
  const [dropDown, setDropDown] = useState(false);

  const onMouseEnter = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(true);
  };
  const onMouseLeave = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(false);
  };

  const { url } = useRouteMatch();

  return (
    <>
      <nav className={sidebarShown ? styles.navbarWithCollapsedSidebar : styles.navbar}>
        <li className={styles.navItem}>
          <NavLink to={`${url}/dashboard`}>
            <img
              alt="gene-icon"
              src={geneIcon}
              style={{ height: '45px', paddingLeft: '15px' }}
            />
          </NavLink>
        </li>

        <ul
          className={styles.navMenu}
          style={{ left: '50%', margin: 'auto' }}
        >
          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to={`${url}/dashboard`}
            >
              Dashboard
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to={`${url}/documentation`}
            >
              Documentation
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to={`${url}/help`}
            >
              Help
            </NavLink>
          </li>
        </ul>

        <li
          className={styles.navItem}
          onMouseEnter={onMouseEnter}
          onMouseLeave={onMouseLeave}
        >
          <NavLink
            className={styles.profileSettings}
            to={`${url}/settings`}
          >
            {`Hi, ${user.firstName}!`}
            <img
              alt="profiledefault"
              src={profiledefault}
              style={{ height: '40px', paddingLeft: '15px' }}
            />
          </NavLink>
          { dropDown && <SettingsDropdown setUser={setUser} /> }
        </li>

      </nav>
    </>

  );
}

export default NavigationBar;
