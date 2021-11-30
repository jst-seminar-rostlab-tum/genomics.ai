import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { IconContext } from 'react-icons/lib';
import * as IoIcons from 'react-icons/io';
import * as CGIcons from 'react-icons/cg';
import styles from './subMenu.module.css';

function SubMenu({ item }) {
  const [subnav, setSubnav] = useState(false);
  const showSubnav = () => setSubnav(!subnav);

  function getDropdownStatusIcon() {
    if (item.subNav && subnav) {
      return (
        <>
          <IoIcons.IoIosArrowDown />
        </>
      );
    }
    return (
      <>
        <IoIcons.IoIosArrowForward />
      </>
    );
  }

  return (
    <IconContext.Provider value={{ color: '#fff' }}>
      <div className={styles.divider} />
      <div
        className={styles.sidebarLink}
        onClick={item.subNav && showSubnav}
      >
        <span className={styles.sidebarLabel}>{item.name}</span>
        <div>
          {getDropdownStatusIcon()}
        </div>
      </div>
      {
            subnav && item.subNav.map((subItem) => (
              <Link
                className={styles.dropdownLink}
                to={subItem.path}
                key={subItem.id}
              >
                <CGIcons.CgFileDocument />
                <span className={styles.sidebarLabel}>
                  {subItem.name}
                </span>
              </Link>
            ))
          }
    </IconContext.Provider>
  );
}

export default SubMenu;
