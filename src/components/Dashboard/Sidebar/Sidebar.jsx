import React, { useState } from 'react';
import { IconContext } from 'react-icons/lib';
import { SidebarData } from './SidebarData';
import SubMenu from './SubMenu/SubMenu';
import styles from './sidebar.module.css';
import arrowOpen from '../../../assets/arrow-open.png';
import arrowClose from '../../../assets/arrow-close.png';
import geneIcon from '../../../assets/gene.png';

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
                              <div
                                className={styles.toggleButton}
                              >
                                <input
                                  type="image"
                                  alt="toggle-icon"
                                  src={arrowOpen}
                                  style={{ height: '30px', transition: '350ms' }}
                                  onClick={showSidebar}
                                />
                              </div>
                            </div>
                          </div>

                          <div style={{ paddingInline: '10px' }}>
                            <input
                              className={styles.projectsBannerCollapsed}
                              type="image"
                              alt="toggle sidebar"
                              src={geneIcon}
                              style={{
                                height: '35px', flex: 1, flexDirection: 'row', alignItems: 'center',
                              }}
                              onClick={showSidebar}
                            />
                          </div>

                        </div>

                      </div>
                    ) : (
                      <div className={styles.sidebarNav}>
                        <div className={styles.sidebarWrap}>
                          <div style={{ flex: 1, flexDirection: 'row', alignItems: 'center' }}>

                            <div style={{ paddingBlock: '25px', paddingLeft: '270px' }}>
                              <input
                                type="image"
                                className={styles.toggleButton}
                                alt="toggle-icon"
                                src={arrowClose}
                                style={{ height: '30px', transition: '350ms' }}
                                onClick={showSidebar}
                              />
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
                          { SidebarData.map((item) => <SubMenu item={item} key={item.id} />)}
                        </div>
                      </div>
                    )
          }
      </IconContext.Provider>
    </>
  );
};

export default Sidebar;
