import React, { useState } from 'react';
import { Link } from 'react-router-dom';
import { IconContext } from 'react-icons/lib';
import { SidebarData } from './SidebarData';
import SubMenu from './SubMenu/SubMenu';
import styles from './sidebar.module.css';
import arrowOpen from '/Users/aminbensaad/GeneCruncher-Front/src/assets/arrow-open.png';
import arrowClose from '/Users/aminbensaad/GeneCruncher-Front/src/assets/arrow-close.png';
import geneIcon from '/Users/aminbensaad/GeneCruncher-Front/src/assets/gene.png';

const Sidebar = () => {
  const [sidebar, setSidebar] = useState(false);
  const showSidebar = () => setSidebar(!sidebar);

  return (
    <>
      <IconContext.Provider value={{ color: '#fff' }}>
        {
                    sidebar ? (
                      <div className={styles.sidebarNavCollapsed}>
                        <div className={styles.sidebarWrap}>
                          <div style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}>

                            <div style={{ paddingBlock: '25px', paddingLeft: '80px' }}>
                              <Link
                                className={styles.toggleButton}
                                to="#"
                              >
                                <img
                                  alt="toggle-icon"
                                  src={arrowOpen}
                                  style={{ height: '30px', transition: '350ms' }}
                                  onClick={showSidebar}
                                />
                              </Link>
                            </div>
                          </div>

                          <div style={{ paddingInline: '10px' }}>
                            <div
                              className={styles.projectsBannerCollapsed}
                              style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}
                              onClick={showSidebar}
                            >
                              <img
                                alt="gene-icon"
                                src={geneIcon}
                                style={{ height: '35px' }}
                              />
                            </div>
                          </div>

                        </div>

                      </div>
                    ) : (
                      <div className={styles.sidebarNav}>
                        <div className={styles.sidebarWrap}>
                          <div style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}>

                            <div style={{ paddingBlock: '25px', paddingLeft: '270px' }}>
                              <Link
                                className={styles.toggleButton}
                                to="#"
                              >
                                <img
                                  alt="toggle-icon"
                                  src={arrowClose}
                                  style={{ height: '30px', transition: '350ms' }}
                                  onClick={showSidebar}
                                />
                              </Link>
                            </div>
                          </div>

                          <div
                            className={styles.projectsBanner}
                            style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}
                          >
                            <img
                              alt="gene-icon"
                              src={geneIcon}
                              style={{ height: '35px', paddingLeft: '5px', paddingRight: '30px' }}
                            />
                            Projects
                          </div>
                          { SidebarData.map((item, index) => <SubMenu item={item} key={index} />)}
                        </div>
                      </div>
                    )
          }
      </IconContext.Provider>
    </>
  );
};

export default Sidebar;
