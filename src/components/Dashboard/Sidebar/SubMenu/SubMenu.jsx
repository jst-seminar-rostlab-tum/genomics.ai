import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { IconContext } from 'react-icons/lib';
import { IoIosArrowDown, IoIosArrowForward, CgFileDocument } from 'react-icons/all';
import styles from './subMenu.module.css';

function SubMenu({ item }) {
  const [subnav, setSubnav] = useState(false);
  const showSubnav = () => setSubnav(!subnav);

  function getDropdownStatusIcon() {
    if (item.subNav && subnav) {
      return (
        <>
          <IoIosArrowDown />
        </>
      );
    }
    return (
      <>
        <IoIosArrowForward />
      </>
    );
  }

  return (
    <IconContext.Provider value={{ color: '#fff' }}>
      <div className={styles.divider} />
      <div
        tabIndex={0}
        onKeyPress={() => {}}
        role="button"
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
            <CgFileDocument />
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
