import React, { useState } from 'react';
import { Link as NavLink } from 'react-router-dom';
import { SettingsDropdown } from './SettingsDropdown/SettingsDropdown';
import styles from './navigationBar.module.css';
import geneIcon from '../../../assets/logo-blue.png';
import profiledefault from '../../../assets/user.png';

function NavigationBar({ sidebarShown }) {
  const [dropDown, setDropDown] = useState(false);

  const onMouseEnter = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(true);
  };
  const onMouseLeave = () => {
    if (window.innerWidth < 960) setDropDown(false);
    else setDropDown(false);
  };

  return (
    <>
      <nav className={sidebarShown ? styles.navbarWithCollapsedSidebar : styles.navbar}>
        <li className={styles.navItem}>
          <NavLink to="/dashboard">
            <img
              alt="gene-icon"
              src={geneIcon}
              style={{ height: '45px', paddingLeft: '15px' }}
              hover={sidebarShown}
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
              to="/dashboard"
            >
              Dashboard
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to="/documentation"
            >
              Documentation
            </NavLink>
          </li>

          <li className={styles.navItem}>
            <NavLink
              className={styles.navLinks}
              to="/help"
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
            to="/settings"
          >
            Amin Ben Saad
            <img
              alt="profiledefault"
              src={profiledefault}
              style={{ height: '40px', paddingLeft: '15px' }}
            />
          </NavLink>
          { dropDown && <SettingsDropdown /> }
        </li>

      </nav>
    </>

  );
}

export default NavigationBar;
